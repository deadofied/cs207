#include <fstream>
#include <time.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <algorithm>
#include <chrono>
#include <thread>
#include <queue>

#include "CS207/MySDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Data.hpp"
#include "Graph.hpp"
#include "Point.hpp"

/** @file Project.cpp
 * @brief Define the functions for various force directed drawing algorithms. 
 */

//=========
//CONSTANTS
//=========

float node_lim_size = 15.0; // Node limit size for SDL Viewer
float edge_lim_size = 1.0; //Edge limit size for SDL Viewer
static constexpr double min_win_size = -0.5;// min x and y coordinates for the random creation of points
static constexpr double max_win_size = 0.5;// max x and y coordinates for the random creation of points

double grav = 4.775;//Gravity (Gauss-Siedel)
double alpha = 0.1; //Temperature (Gauss-Siedel) 
double charge = -3.0;//Charge (Gauss-Siedel)
double linkDistance = 1.0;//Minimum Edge distance (Gauss-Siedel)
double linkStrength = 1.0;// Edge strength (Gauss-Siedel)
double friction = 0.9;//Friction(Gauss-Siedel)
double theta2 = 0.64;//Distance threshold (Gauss-Siedel)

float t = 0.9; // Temperature of the system (Fruchterman-Reingold)

//Define the Center of Gravity of the SDL Viewer (Gauss-Siedel)
Point center = Point((double)min_win_size + (double)(max_win_size-min_win_size)/2.0,(double)min_win_size + (double)(max_win_size-min_win_size)/2.0,0.0);

//=======================
//FUNCTORS FOR SDL VIEWER
//=======================


/** Node size function object for use in the SDLViewer. 
 * This functor draw sizes that are proportional to the degree of a given node
*/
struct NodeSize {
 template <typename NODE>
 float operator()(const NODE& n) {
 
   return n.value().mass;//3*n.degree();

 }
};

/** Node size function object for use in the SDLViewer (Fruchterman-Reingold)
 * This functor draw sizes that are proportional to the degree of a given node
 */
struct NodeSizeStd {
 template <typename NODE>
 float operator()(const NODE& n) {
   (void)n;
   return 10.0;
 }
};

/** Edge Width function object for use in the SDLViewer (Gauss-Siedel) */
struct EdgeWidth {
  EdgeWidth(){}

 template <typename EDGE>
 float operator()(const EDGE& e) {
 
   if (e.value()!=0.0) return e.value()/4.0;
   else return 1.0;

 }
};

/** Edge Width function object for use in the SDLViewer (Fruchterman Reingold) */
struct EdgeWidthStd {
 template <typename EDGE>
 float operator()(const EDGE& e) {
   (void)e;
   return 0.5;

 }
};

/**Color Functor to pass to the SDL Viewer.*/
struct MyColor {
 template <typename NODE>
 CS207::Color operator()(const NODE& node) {
   if (node.value().group == 0) return CS207::Color(102.0/255.0,51.0/255.0,0.0/255.0);//Brown
   else if (node.value().group == 1) return CS207::Color(51.0/255.0,51.0/255.0,255.0/255.0);//Blue
   else if (node.value().group == 2) return CS207::Color(255.0/255.0,0.0/255.0,0.0/255.0);//Red
   else if (node.value().group == 3) return CS207::Color(0.0/255.0,102.0/255.0,0.0/255.0);//Green
   else if (node.value().group == 4) return CS207::Color(255.0/255.0,128.0/255.0,0.0/255.0);//Orange
   else if (node.value().group == 5) return CS207::Color(153.0/255.0,51.0/255.0,255.0/255.0);//Violet
   else if (node.value().group == 6) return CS207::Color(0.0/255.0,204.0/255.0,204.0/255.0);//Sky Blue
   else if (node.value().group == 7) return CS207::Color(255.0/255.0,178.0/255.0,102.0/255.0);//Skin
   else if (node.value().group == 8) return CS207::Color(178.0/255.0,255.0/255.0,102.0/255.0);//Bright green
   else if (node.value().group == 9) return CS207::Color(255.0/255.0,153.0/255.0,204.0/255.0);//Pink
   else return CS207::Color(160.0/255.0,160.0/255.0,160.0/255.0);//Grey
 }
};

/**Color Functor to pass to the SDL Viewer (Fruchterman Reingold).*/
struct MyColorFruchterman {
 template <typename NODE>
 CS207::Color operator()(const NODE& node) {
   (void)node;
   return CS207::Color::make_heat(float(node.index())/5.0);
 }
};

/**Color Functor to pass to the SDL Viewer.*/
struct MyColorKKSA {
 template <typename NODE>
 CS207::Color operator()(const NODE& node) {
   if      (node.value().color == 1) return CS207::Color(51.0/255.0,51.0/255.0,255.0/255.0);    //Blue
   else if (node.value().color == 2) return CS207::Color(255.0/255.0,0.0/255.0,0.0/255.0);      //Red
   else if (node.value().color == 8) return CS207::Color(178.0/255.0,255.0/255.0,102.0/255.0);  //Bright green
   else return CS207::Color::make_heat(1);
 }
};

//========================================
//SHORTEST PATH USING BREADTH-FIRST SEARCH
//========================================

template <typename NODE, typename GRAPH>
float graph_dist(NODE n1, NODE n2, GRAPH& g) {
  for (auto n = g.RA_node_begin(); n != g.RA_node_end(); ++n) { (*n).value().graphdist = -1; }

  (g.node(n1.index())).value().graphdist = 0;

  std::queue<NODE> node_q;
  node_q.push(g.node(n1.index()));

  // search graph and update values
  while (!node_q.empty()) {
    auto i = node_q.front();  // remove from queue
    for (auto adj = i.edge_begin(); adj != i.edge_end(); ++adj) { // for each adjacent edge

      NODE adjn;

      if ((*adj).node1() == i)
        adjn = (*adj).node2();
      else
        adjn = (*adj).node1();

      if (adjn.value().graphdist == -1) {       // if never been reached

        adjn.value().graphdist = i.value().graphdist + 1; // add value
        node_q.push(adjn);             // put this adj node into queue
      }
    }

    node_q.pop();
  }

  return (g.node(n2.index())).value().graphdist;
}


//==============
//STEP FUNCTIONS
//==============


/** Perform a step in the Gauss siedel algorithm that brings us one step closer to convergence
 *
 * @param[in,out]   g      valid Graph object 
 * @param[in]  vector of node positions
 *
 */
template <typename G,typename V>
void Gauss_Siedel(G& g,V px_s) {

    //STEP 1: GAUSS-RIEDEL RELAXATION FOR LINKS
    for(auto ei=g.edge_begin();ei!=g.edge_end();++ei){
      auto edge = *ei;
      auto s = edge.node1();
      auto t = edge.node2();
      Point vec = t.position() - s.position();
      double distance = alpha * linkStrength * (norm(vec)-linkDistance) / pow(norm(vec),2.0);
      vec*= distance;
      double k = ((double)s.degree())/((double)s.degree()+(double)t.degree());
      t.set_position(t.position() - vec * k);
      s.set_position(s.position() - vec * (1.0-k));
    }

    //STEP 2: GRAVITY FORCES
    for (auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      n.set_position( n.position() + (center - n.position()) * alpha * grav);
    }

    //STEP 3: COMPUTE QUADTREE CENTER OF MASS AND APPLY CHARGE FORCES
    for(auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      for(auto ni2=g.RA_node_begin();ni2!=g.RA_node_end();++ni2){
        auto n2 = *ni2;
        double dx = n.position().x - n2.position().x;
        double dy = n.position().y - n2.position().y;
        double dn = pow(dx,2.0) + pow(dy,2.0);
        double dw = max_win_size - min_win_size;

        if(pow(dw,2.0)/theta2<dn){
          px_s[n2.index()].x -= dx * charge/dn;
          px_s[n2.index()].y -= dy * charge/dn;
        }
      }

    }

    //STEP 4: POSITION VERLET INTEGRATION
    for (auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      n.set_position( n.position() - (px_s[n.index()] - n.position()) * friction);
      px_s[n.index()] = n.position();
    }
}

/** Perform a step in the Gauss siedel algorithm that brings us one step closer to convergence
 *
 * @param[in,out]   g      valid Graph object 
 * @param[in]  vector of node positions
 *
 */
template <typename G,typename V>
void Gauss_Siedel_with_Subgroups(G& g,V px_s) {

    grav = 5.5;
    //STEP 1: GAUSS-RIEDEL RELAXATION FOR LINKS
    for(auto ei=g.edge_begin();ei!=g.edge_end();++ei){
      auto edge = *ei;
      auto s = edge.node1();
      auto t = edge.node2();
      Point vec = t.position() - s.position();
      double distance = alpha * linkStrength * (norm(vec)-linkDistance) / pow(norm(vec),2.0);
      vec*= distance;
      double k = ((double)s.degree())/((double)s.degree()+(double)t.degree());
      t.set_position(t.position() - vec * k);
      s.set_position(s.position() - vec * (1.0-k));
    }

    //STEP 2: GRAVITY FORCES
    for (auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      n.set_position( n.position() + (center - n.position()) * alpha * grav);
    }

    //STEP 3: COMPUTE QUADTREE CENTER OF MASS AND APPLY CHARGE FORCES
    for(auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      for(auto ni2=g.RA_node_begin();ni2!=g.RA_node_end();++ni2){
        auto n2 = *ni2;
        double dx = n.position().x - n2.position().x;
        double dy = n.position().y - n2.position().y;
        double dn = pow(dx,2.0) + pow(dy,2.0);
        double dw = max_win_size - min_win_size;

        if(pow(dw,2.0)/theta2<dn && n.value().group!=n2.value().group){
          px_s[n2.index()].x -= dx * charge/dn;
          px_s[n2.index()].y -= dy * charge/dn;
        }
      }

    }

    //STEP 4: POSITION VERLET INTEGRATION
    for (auto ni=g.RA_node_begin();ni!=g.RA_node_end();++ni){
      auto n = *ni;
      n.set_position( n.position() - (px_s[n.index()] - n.position()) * friction);
      px_s[n.index()] = n.position();
    }
}

 /** Perform a step in the Fruchterman Reingold algorithm that brings us one step closer to convergence
 *
 * @param[in,out]   g      valid Graph object 
 * @param[in]     t     temperature of system at the current state
 * @pre  0 < t <= 1
 * @pre  g.num_nodes() > 2
 * In this algorithm, the nodes are represented by steel rings and the edges are springs between them. The attractive force is 
 * analogous to the spring force and the repulsive force is analogous to the electrical force. The basic idea is to minimize the 
 * energy of the system by moving the nodes and changing the forces between them. For more details refer to the Force Directed algorithm.
 * In this algorithm, the sum of the force vectors determines which direction a node should move. The step width, which is a constant 
 * determines how far a node moves in a single step. When the energy of the system is minimized, the nodes stop moving and the system 
 * reaches it's equilibrium state. The drawback of this is that if we define a constant step width, there is no guarantee that the system 
 * will reach equilibrium at all. T.M.J. Fruchterman and E.M. Reingold introduced a "global temperature" that controls the step width of 
 * node movements and the algorithm's termination. The step width is proportional to the temperature, so if the temperature is hot, the nodes 
 * move faster (i.e, a larger distance in each single step). This temperature is the same for all nodes, and cools down at each iteration. 
 * Once the nodes stop moving, the system terminates.
 */
template <typename G>
void Fruchterman_Reingold(G& g, float& t){
  float k = sqrt(0.5/g.num_nodes()); //Normalizing constant (Fruchterman-Reingold)

    //Sleep Time
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // Repulsive Forces on electrons
    for (auto node1 = g.RA_node_begin(); node1 != g.RA_node_end(); ++node1) { 
      auto n1 = *node1;

      for (auto node2 = g.RA_node_begin(); node2 != g.RA_node_end(); ++node2) {
        auto n2 = *node2;

        if (n1 != n2) {
          Point diff = n1.position() - n2.position();
          n1.value().displacement = n1.value().displacement + diff / norm(diff) * ((k * k) / norm(diff));
        } 
      }
    }

    // Calculate the attractive forces on edges/springs
    for (auto e = g.edge_begin(); e != g.edge_end(); ++e) {
      auto edge = *e;

      Point diff = edge.node1().position() - edge.node2().position();
      edge.node1().value().displacement -= ((diff / norm(diff)) * (norm(diff) * norm(diff) / k));
      edge.node2().value().displacement += ((diff / norm(diff)) * (norm(diff) * norm(diff) / k));

    }

    for (auto node = g.RA_node_begin(); node != g.RA_node_end(); ++node) {
      auto n = *node;

      Point pos = n.position() + (n.value().displacement/norm(n.value().displacement)) * std::min(norm(n.value().displacement), (double)t);
      
      pos.x = std::min(0.5, std::max(-0.5, pos.x));
      pos.y = std::min(0.5, std::max(-0.5, pos.y)); 
      pos.y = std::min(0.5, std::max(-0.5, pos.y));

      n.set_position(pos);

      t = .9 * t;
    //}
  }
}

/* Updates the step size, necessary for the simulated annealing Kamada Kawai approach 
 * @post: @param progress has been updated
 */
Point update_step_length(Point step, float old_E, float new_E, float t, float& progress) {
  if (new_E < old_E) {
    progress += 1;
    if (progress >= 5) {
      progress = 0;
      step = step / t;
    }
  }
  else {
    progress = 0;
    step = t * step;
  }
  return step;
}

/* The Kamada Kawai simulated annealing algorithm is best suited for social networks.
 * It gradually gets closer to a "good" representation of a network (one in which the
 * euclidean distance is as close as possible to the graph distance. Because Kamada &
 * Kawai's algorithm is particularly expensive - O(num_nodes^4) - an approach close to
 * simulated annealing was used here.
 * @pre: 0 < t < 1
 * @pre: step != Point(0,0,0)
 * @post:the nodes of @param g have moved a step closer to the Kamada Kawai optimal distribution
 */
template <typename G>
void KamadaSA(G& g, float& energy, float& energy0, float& progress, float& t, Point& step){
    std::this_thread::sleep_for(std::chrono::milliseconds(250));

    energy0 = energy;
    energy = 0;
    for (auto i = g.RA_node_begin(); i != g.RA_node_end(); ++i) {
      Point f = Point(0,0,0);

      //update the force for adjacent nodes
      auto nod1 = (*i);
      for (auto edg = nod1.edge_begin(); edg != nod1.edge_end(); ++edg) {
        float dist = (*edg).length();
        f += ( dist - graph_dist( (*i),(*edg).node2(),g ) ) / dist * (((*edg).node2()).position() - (*i).position());
      }

      //updates the force for all nodes
      for (auto j = g.RA_node_begin(); j != g.RA_node_end(); ++j) {
        if (i != j) {
          float dist = pow( pow((*i).position().x - (*j).position().x, 2) +
                            pow((*i).position().y - (*j).position().y, 2) +
                            pow((*i).position().z - (*j).position().z, 2), .5);
          f += ( dist - graph_dist( (*i),(*j),g ) ) / dist * ((*j).position() - (*i).position());
        }
      }
      (*i).position() += step * (f / norm(f));
      energy += normSq(f);
    } 
    step = update_step_length(step, energy0, energy, t, progress);
}


//=============
//MAIN FUNCTION
//=============


int main() {

    int algo;
    std::cout << "Please select a network force algorithm:" << std::endl;
    std::cout << "(1) Gauss-Seidel, stabilize in concentric circles (Les Miserables dataset) " << std::endl;
    std::cout << "(2) Gauss-Seidel, stabilize in groups (Les Miserables dataset)" << std::endl;
    std::cout << "(3) Fruchterman-Reingold " << std::endl;
    std::cout << "(4) Kamada Kawai with Simulated Annealing Approach " << std::endl;
    std::cout << "(5) Kamada Kawai (Les Miserables dataset) " << std::endl;
    std::cout << "(6) Gauss-Seidel, stabilize in groups (US Senators voting dataset)" << std::endl;
    std::cout << "(7) Kamada Kawai (US Senators voting dataset)" << std::endl;
    std::cout << "(8) Kamada Kawai (karate dataset)" << std::endl;
    std::cin >> algo;

    do {
      if (algo != 1 && algo != 2 && algo != 3 && algo != 4 && algo != 5 && algo != 6 && algo != 7 && algo != 8) {
        std::cout << "Invalid selection. Please try again." << std::endl;
        std::cin >> algo;
      }
    } while (!(algo == 1 || algo == 2 || algo == 3 || algo == 4 || algo == 5 || algo == 6 || algo == 7 || algo == 8));


//========
//TYPEDEFS
//========


// Define your Graph type; this is a placeholder.
typedef Graph<NodeValue,float> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;

typedef Graph<NodeValueFR,float> GraphFR;
typedef GraphFR::Node NodeFR;
typedef GraphFR::Edge EdgeFR;

typedef Graph<NodeValueKKSA,float> GraphKKSA;
typedef GraphKKSA::Node NodeKKSA;
typedef GraphKKSA::Edge EdgeKKSA;

  GraphType g1;
  GraphFR g2;
  GraphKKSA g3;

  //Gauss
  double m = 15.0;

  //Kamada SA
  float t = 0.9;
  Point step = Point(5,5,5);
  float energy = std::numeric_limits<float>::max();
  float energy0;
  float progress = 0;

  if (algo == 1 || algo == 2) initialize_les_miserables(g1,min_win_size,max_win_size,m);
  else if (algo == 3) initialize_fruchterman(g2);
  else if (algo == 4) initialize_square(g3);
  else if (algo == 5) initialize_les_miserables_KKSA(g3);  
  else if (algo == 6) initialize_senate(g1,min_win_size,max_win_size,m);
  else if (algo == 7) initialize_senate_kk(g3);
  else if (algo == 8) initialize_karate(g3);

  // Print out the statistics
  //std::cout << g.num_nodes() << " " << g.num_edges() << std::endl;
  
  //Initialize px_s (Gauss_Siedel)
  std::vector<Point> px_s;
  for (auto ni = g1.node_begin(); ni!=g1.node_end(); ++ni){
    px_s.push_back((*ni).position());
  }


  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map1 = viewer.empty_node_map(g1);
  auto node_map2 = viewer.empty_node_map(g2);
  auto node_map3 = viewer.empty_node_map(g3);

  viewer.launch();

  if (algo == 1 || algo == 2 || algo == 6) {
    viewer.add_nodes(g1.node_begin(), g1.node_end(),MyColor(),CS207::DefaultPosition(), NodeSize(), node_map1);
    viewer.add_edges(g1.edge_begin(), g1.edge_end(), EdgeWidth(), node_map1);
  }
  else if (algo == 3) {
    viewer.add_nodes(g2.node_begin(), g2.node_end(),MyColorFruchterman(),CS207::DefaultPosition(), NodeSizeStd(), node_map2);
    viewer.add_edges(g2.edge_begin(), g2.edge_end(), EdgeWidthStd(), node_map2);
  }
  else if (algo == 4 || algo == 5 || algo == 7 || algo == 8) {
    viewer.add_nodes(g3.node_begin(), g3.node_end(), MyColorKKSA(), CS207::DefaultPosition(), NodeSizeStd(), node_map3);
    viewer.add_edges(g3.edge_begin(), g3.edge_end(), EdgeWidthStd(), node_map3);
  }
  
  viewer.center_view();

  // Begin the simulation
  while(true) {
    // Call function to perform a step of a force-directed algorithm
    if (algo == 1) Gauss_Siedel(g1,px_s);
    else if (algo == 2 || algo == 6) Gauss_Siedel_with_Subgroups(g1,px_s);
    else if (algo == 3) Fruchterman_Reingold(g2, t);
    else if (algo == 4 || algo == 5 || algo == 7 || algo == 8) KamadaSA(g3, energy, energy0, progress, t, step);
    
    // Update viewer with nodes' new positions
    if (algo == 1 || algo == 2 || algo == 6) viewer.add_nodes(g1.node_begin(), g1.node_end(),MyColor(),CS207::DefaultPosition(), NodeSize(), node_map1);
    else if (algo == 3) viewer.add_nodes(g2.node_begin(), g2.node_end(), MyColorFruchterman(),CS207::DefaultPosition(), NodeSizeStd(), node_map2);
    else if (algo == 4 || algo == 5 || algo == 7 || algo == 8) {
      viewer.add_nodes(g3.node_begin(), g3.node_end(), MyColorKKSA(), CS207::DefaultPosition(), NodeSizeStd(), node_map3);
      viewer.center_view();
    }

    viewer.set_label(t);

    if (algo == 1 || algo == 2)
      if (g1.size() < 100)
        CS207::sleep(0.01);
  }
  return 0;
}
