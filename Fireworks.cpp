/** owner of this c++ file is nikhil :) **/


#include <fstream>
#include <time.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <algorithm>


#include "CS207/MySDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"

#include <chrono>
#include <thread>

using namespace std;

//===============
//NODE VALUE TYPE
//===============

struct NodeValue {
  NodeValue() : displacement(Point(0.0, 0.0, 0.0)), neighbors(0) {}
  Point displacement;
  int neighbors;
};

//========
//TYPEDEFS
//========

// Define your Graph type; this is a placeholder.
typedef Graph<NodeValue,float> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;


unsigned n_nodes = 500;

//=======================
//FUNCTORS FOR SDL VIEWER
//=======================



/** Node size function object for use in the SDLViewer. 
 * This functor draw sizes that are proportional to the degree of a given node
*/
struct NodeSize {
 template <typename NODE>
 int operator()(const NODE& n) {
  (void)n;
   return 10.0;
 }
};

/** Edge Width function object for use in the SDLViewer. 
*/
struct EdgeWidth {
 template <typename EDGE>
 int operator()(const EDGE& e) {
 
   if (e.value()!=0.0) return e.value();
   else return 0.1;

 }
};

/**Color Functor to pass to the SDL Viewer.*/
struct MyColor {
 template <typename NODE>
 CS207::Color operator()(NODE& node) {
   return CS207::Color::make_heat((float(float(node.value().neighbors))/n_nodes));
 }
};


//=============
//MAIN FUNCTION
//=============


int main() {

  GraphType g;
  srand (time(NULL));
  static std::mt19937 default_generator;
  std::uniform_real_distribution<double> dist(0.0, 1.0);



std::vector<Node> nod_;
  for(unsigned i=0; i<n_nodes;++i){
    //float m = node_lim_size*dist(default_generator);
    double x = dist(default_generator);
    double y = dist(default_generator);
    double z = dist(default_generator);
    //NodeValue nv = NodeValue(m);
    nod_.push_back(g.add_node(Point(x,y,z),NodeValue()));
  }





  // Print out the statistics
  std::cout << g.num_nodes() << " " << g.num_edges() << std::endl;

  

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(g);
  viewer.launch();

  //viewer.add_nodes(g.node_begin(), g.node_end(), node_map);
  viewer.add_nodes(g.node_begin(), g.node_end(), MyColor(), CS207::DefaultPosition(), NodeSize(), node_map);
  //viewer.add_edges(g.edge_begin(), g.edge_end(), EdgeWidth(), node_map);
  viewer.center_view();




  for (int i = 0; i < 5000; ++i) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(100));


    // gravity
    for (auto n = g.node_begin(); n != g.node_end(); ++n) {
      auto node = *n;

      node.value().displacement -= 0.001 * node.position();
      node.set_position(node.position() + node.value().displacement);
    } 


    // calc neighbors, used for colors
    for (auto n1 = g.node_begin(); n1 != g.node_end(); ++n1) {

      auto node1 = *n1;
      node1.value().neighbors = 0;

      for (auto n2 = g.node_begin(); n2 != g.node_end(); ++n2) {
        
        auto node2 = *n2;

        if (node1 != node2) {
          if (norm(node2.position() - node1.position()) < 0.25)
            node1.value().neighbors+=1;
        }

      }
    } 

    viewer.add_nodes(g.node_begin(), g.node_end(), MyColor(), CS207::DefaultPosition(), NodeSize(), node_map);
    viewer.set_label(i);
  }

  return 0;
}
