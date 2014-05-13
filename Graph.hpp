#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"
//#include "Mesh.hpp"

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> 
class Graph {
 public:

  // PUBLIC TYPE DEFINITIONS

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;

  typedef unsigned node_id_type;

  typedef V node_value_type;

  typedef E edge_value_type;

 private:
  struct internal_node {
   Point position_;
   node_value_type value_;
   size_type index_; 

  internal_node(Point position, node_value_type value, size_type index) : position_(position), value_(value), index_(index) {
   }
  };

  struct internal_edge {
    node_id_type uid_;  
    edge_value_type value_; 

    // define some operators so that we can iterate through edges nicely
    bool operator==(const internal_edge& e) const {
      return (uid_ == e.uid_);
    }
    bool operator<(const internal_edge& e) const {
        return (uid_ < e.uid_);
    }
    bool operator>(const internal_edge& e) const {
        return (uid_ > e.uid_);
    }
    bool operator!=(const internal_edge& e) const {
        return (uid_ != e.uid_);
    }
    bool operator<=(const internal_edge& e) const {
        return (uid_ <= e.uid_);
    }
    bool operator>=(const internal_edge& e) const {
        return (uid_ >= e.uid_);
    }

    internal_edge(node_id_type uid, edge_value_type value = edge_value_type()) : uid_(uid), value_(value) {
    }
  };

 public:
  /** Type of node iterators, which iterate over all graph nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class incident_iterator;

  // CONSTRUCTOR AND DESTRUCTOR

  Graph() : nodes_(), i2u_(), bridge_() {

  }

  /** Default destructor */

  ~Graph() = default;

  // NODES

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
  //class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {

    }

    /** Return this node's position. */
    Point& position() const {
      return graph_->nodes_[uid_].position_;
    }

    /** Set this node's position. */
    void set_position ( const Point & p) {
      graph_->nodes_[uid_].position_ = p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].index_;
    }

    /** Return the current node's value. */
    node_value_type& value() {
      return graph_->nodes_[uid_].value_;
    }

    /** Same as above, can't modify the output. */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value_;      
    }

    /** Test for equality between two nodes.
     * @param[in] @n node of a graph
     * @pre @a n is a valid node of a graph
     * @return true if some @a n has the same index as the current node and if a is a valid node of this graph
     */
    bool operator==(const Node& n) const {
        return (uid_ == n.uid_ &&  graph_ == n.graph_);
    }

    /** Test whether one node is smaller than another
     * @param[in] @n node of this graph
     * @pre @a n is a valid node of this graph
     * @return true if some @a n has a smaller index than the current node 
     */
    bool operator< (const Node& n) const {
        return (uid_ < n.uid_);
    }


    /** Return the degree of the current node */
    size_type degree() const {
      return graph_->bridge_[uid_].size();
    }

    /** Return an incident_iterator for the current node.
    * @post the incident_iterator is initialized
    */
    incident_iterator edge_begin() const {
      // this is the ID of my first edge
      return incident_iterator(graph_, uid_, graph_->bridge_[uid_].begin());
    }
    
    /** Return an incident_iterator for the current node.
    * This function defines the end for the incident_iterator
    */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, uid_, graph_->bridge_[uid_].end());

    }

   private:
    // Only Graph can access our private members
    friend class Graph;
    graph_type* graph_;
    node_id_type uid_;
    
    Node(const graph_type* graph, node_id_type uid): graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size(); // number of active nodes
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */

  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_node(position, value, i2u_.size());
    nodes_.push_back(new_node);
    
    i2u_.push_back(nodes_.size() - 1);
    
    // create empty set at the end of the bridge to store the other nodes that connect to this node
    std::set<internal_edge> temp_set;
    bridge_.push_back(temp_set);

    return Node(this, nodes_.size() - 1);   
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= i2u_.size())
      return false;
    else
      return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < size()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) { 
    assert(i < i2u_.size());
    return Node(this, i2u_[i]);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph.
   * @post new size() == old size() - 1
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Polynomial in size().
   */
  void remove_node(const Node& n) {
    size_type index = n.index();
    node_id_type uid = i2u_[index];

    // empty the set which shows who this node is connected to
    bridge_[uid].clear();

    // remove all connections to this node from the bridge, i.e. break all the edges in the bridge
    internal_edge new_edge(n.uid_);    

    for(node_id_type i = 0; i < bridge_.size(); i++) {
      bridge_[i].erase(new_edge);
    }

    // remove the node
    i2u_.erase(i2u_.begin() + index);

    // now need to update all the indices of the nodes
    for(size_type i = 0; i < i2u_.size(); i++)
        nodes_[i2u_[i]].index_ = i;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    bridge_.clear();
  }


  // EDGES

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    Edge(const graph_type* graph, node_id_type node1uid, node_id_type node2uid) : graph_(const_cast<graph_type*>(graph)), node1uid_(node1uid), node2uid_(node2uid) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(graph_->nodes_[node1uid_].index_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(graph_->nodes_[node2uid_].index_);
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return the current edge's value. */
    edge_value_type& value() {
      internal_edge temp_edge(node2uid_);

      auto toReturn = edge_value_type();

      for (auto it = graph_->bridge_[node1uid_].begin(); it != graph_->bridge_[node1uid_].end(); ++it)
          if (*it == temp_edge) toReturn = (*it).value_;

      return toReturn;
    }

    /** Same as above, can't modify the output. */
    const edge_value_type& value() const {
      internal_edge temp_edge(node2uid_);

      for (auto it = graph_->bridge_[node1uid_].begin(); it != graph_->bridge_[node1uid_].end(); ++it)
          if (*it == temp_edge) return (*it).value_;
      assert(0);
      return (*(graph_->bridge_[node1uid_].begin())).value_;
    }

    bool operator==(const Edge& e) const {
        return (((node1() == e.node1() && node2() == e.node2()) || (node1() == e.node2() && node1() == e.node2())) && (graph_ == e.graph_));
    }
    bool operator< (const Edge& e) const {
        return node1uid_ < e.node2uid_;
    }

   private:
    // Only Graph can access our private members
    friend class Graph;
    graph_type* graph_;
    node_id_type node1uid_;
    node_id_type node2uid_;

    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    unsigned count = 0;

    for(node_id_type i = 0; i < bridge_.size(); i++)
      count += bridge_[i].size();
    
    return (count / 2);
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());

    auto ei = edge_begin();
    for (size_type n = 0; n < i ; n++) ++ei;
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    internal_edge temp_edge(b.uid_);

    for (auto it = bridge_[a.uid_].begin() ; it != bridge_[a.uid_].end(); ++it)
        if (*it == temp_edge) return true;

    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {

    internal_edge new_edge1(b.uid_, value);
    internal_edge new_edge2(a.uid_, value);

    bridge_[a.uid_].insert(new_edge1);
    bridge_[b.uid_].insert(new_edge2);
    
    return Edge(this, a.uid_, b.uid_); 
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] a,b The nodes potentially defining an edge to be removed.
   * @return 1 if old has_edge(@a a, @a b), 0 otherwise
   * @pre @a a and @a b are valid nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - result
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node& a, const Node& b) {
    internal_edge new_edge1(b.uid_);    
    internal_edge new_edge2(a.uid_);

    if (bridge_[a.uid_].erase(new_edge1) && bridge_[b.uid_].erase(new_edge2)) return 1;
    else return 0;
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] e The edge to remove
   * @pre @a e is a valid edge of this graph
   * @pre has_edge(@a e.node1(), @a e.node2())
   * @post !has_edge(@a e.node1(), @a e.node2())
   * @post new num_edges() == old num_edges() - 1
   *
   * This is a synonym for remove_edge(@a e.node1(), @a e.node2()), but its
   * implementation can assume that @a e is definitely an edge of the graph.
   * This might allow a faster implementation.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  // ITERATORS

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class node_iterator: totally_ordered<node_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;



    /** Construct an invalid node_iterator. */
    node_iterator(){
    }
    
    /** Return the node object that corresponds to the node iterator's current position */
    Node operator*() const{
      return graph_->node(index_);
    }
    
    /** Increment the index of the node iterator */
    node_iterator& operator++(){
      ++index_;
      return *this;
    }

    size_type get_index() {
      return index_;
    }
       
    /** Test whether two node iterators are equal.
     * @param[in] @ni node_iterator of a graph
     * @pre @a ni is a valid node_iterator of a graph
     * @return true if some @a ni has the same index as the current node_iterator and if ni is a valid node_iterator of this graph
     */       
    bool operator==(const node_iterator& n) const{
      return (index_ == n.index_ && graph_ == n.graph_);
    } 

   private:
    friend class Graph;
    friend class edge_iterator;
    graph_type* graph_;
    size_type index_;

    node_iterator(const graph_type* graph, size_type index): graph_(const_cast<graph_type*>(graph)), index_(index) {
    }

  };

  /** Return a node_iterator for this graph.
    *@post give the starting point of the node_iterator.
    */
  node_iterator node_begin() {
    return node_iterator(this, 0);
  }


  node_iterator node_begin(int x) {
    return node_iterator(this, x);
  }
  
  /** Return a node_iterator for this graph.
    *@post give the end point of the node_iterator.
    */
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }



  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class RA_node_iterator: totally_ordered<node_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;



    /** Construct an invalid node_iterator. */
    RA_node_iterator(){
    }
    
    /** Return the node object that corresponds to the node iterator's current position */
    Node operator*() const{
      return graph_->node(index_);
    }
    
    /** Increment the index of the node iterator */
    RA_node_iterator& operator++(){
      ++index_;
      return *this;
    }

    /** Increment the index of the node iterator */
    RA_node_iterator& operator--(){
      --index_;
      return *this;
    }

    size_type get_index() {
      return index_;
    }
       
    /** Test whether two node iterators are equal.
     * @param[in] @ni node_iterator of a graph
     * @pre @a ni is a valid node_iterator of a graph
     * @return true if some @a ni has the same index as the current node_iterator and if ni is a valid node_iterator of this graph
     */       
    bool operator==(const RA_node_iterator& n) const{
      return (index_ == n.index_ && graph_ == n.graph_);
    }

    bool operator!=(const RA_node_iterator& n) const{
      return !(index_ == n.index_ && graph_ == n.graph_);
    } 

    RA_node_iterator& operator+=(int n) const{
      index_ += n;
      return *this;
    }

    RA_node_iterator& operator-=(int n) const{
      index_ -= n;
      return *this;
    }

    RA_node_iterator& operator[](int n) const{
      index_ += n;
      return *this;
    }

    RA_node_iterator operator+(int n) const{
      return RA_node_iterator(graph_, index_ + n);
    }

    RA_node_iterator operator-(int n) const{
      return RA_node_iterator(graph_, index_ - n);
    }

    int operator-(RA_node_iterator ra_it) const{
      return index_ - ra_it.index_;
    }

    bool operator<(RA_node_iterator rai2) {
      return index_ - rai2.index_ > 0;
    }

    bool operator>(RA_node_iterator rai2) {
      return index_ - rai2.index_ < 0;
    }

    bool operator<=(RA_node_iterator rai2) {
      return !(index_ - rai2.index_ < 0);
    }

    bool operator>=(RA_node_iterator rai2) {
      return !(index_ - rai2.index_ > 0);
    }

   private:
    friend class Graph;
    friend class edge_iterator;
    graph_type* graph_;
    size_type index_;

    RA_node_iterator(const graph_type* graph, size_type index):
            graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  }; 


  RA_node_iterator RA_node_begin() {
    return RA_node_iterator(this, 0);
  }
  
  /** Return a node_iterator for this graph.
    *@post give the end point of the node_iterator.
    */
  RA_node_iterator RA_node_end() const {
    return RA_node_iterator(this, num_nodes());
  }

  /** Return a node_iterator for this graph.
    *@post give the end point of the node_iterator.
    */
  RA_node_iterator RA_node_min() const {
    return RA_node_iterator(this, -1);
  }








  /** @class Graph::edge_iterator
   * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator: totally_ordered<edge_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid edge_iterator. */
    edge_iterator() {
    }

    /** Return the Edge corresponding to this edge_iterator.*/
    Edge operator*() const{
      auto pos = graph_->bridge_[i_].begin();
      std::advance(pos, j_);
      return Edge(graph_, i_, (*pos).uid_);
    }

    /**Increment the edge iterator.*/
    edge_iterator& operator++() {

      ++j_;

      while (j_==graph_->bridge_[i_].size() && i_<graph_->bridge_.size()-1) {
        ++i_;
        
        while(graph_->bridge_[i_].size()==0 && i_<graph_->bridge_.size()-1)
          ++i_;
        j_=0;
        
        auto it = (graph_->bridge_[i_]).begin();
        
        while ((*it)<i_ && it!=(graph_->bridge_[i_]).end()) {
          std::advance(it, 1);
          ++j_;
        }
      }
      return *this;
    }

    /** Test whether two edge iterators are equal.
     * @param[in] @ei edge_iterator of a graph
     * @pre @a ei is a valid edge_iterator of a graph
     * @return true if some @a ei has the same two indices as the current edge_iterator and if ei is a valid edge_iterator of this graph
     */
    bool operator==(const edge_iterator& n) const{
      return (i_ == n.i_ && j_ == n.j_ && graph_ == n.graph_);
    } 

   private:
    friend class Graph;
    graph_type* graph_;
    node_id_type i_;
    node_id_type j_;

    edge_iterator(const graph_type* graph, node_id_type index1, node_id_type index2): graph_(const_cast<graph_type*>(graph)), i_(index1), j_(index2) {
    }  
  };

/** Return an edge_iterator for this graph.
    *@post give the starting point of the edge_iterator.
    */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0, 0);
  }
  
  /** Return an edge_iterator for this graph.
    *@post give the end point of the edge_iterator.
    */
  edge_iterator edge_end() const {
    return edge_iterator(this, bridge_.size()-1, bridge_[bridge_.size()-1].size());
  }


  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator: totally_ordered<incident_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_iterator() {
    }

    /** Return the Edge corresponding to this incident_iterator.*/
    Edge operator*() const {
      assert(nodeuid_ < graph_->bridge_.size());
      return Edge(graph_, nodeuid_, (*position_).uid_);
    }

    /** Increment the incident_iterator.*/
    incident_iterator& operator++() {
      ++position_;
      return *this;
    }

    /** Test whether two incident iterators are equal
     * @param[in] @n incident_iterator of a graph
     * @pre @a n is a valid incident_iterator of a graph
     * @return true if some @a n has the same index and refers to the same node as the current incident_iterator and if n is a valid incident_iterator of this graph
     */
    bool operator==(const incident_iterator& n) const{
      return (position_ == n.position_ && nodeuid_ == n.nodeuid_ && graph_ == n.graph_);  
    } 
    
   private:
    friend class Graph;
    graph_type* graph_;
    node_id_type nodeuid_;
    typename std::set<internal_edge>::iterator position_;
    
    incident_iterator(const graph_type* graph, node_id_type nodeuid, typename std::set<internal_edge>::iterator itr): graph_(const_cast<graph_type*>(graph)), nodeuid_(nodeuid), position_(itr) {}


  };

 private:
  std::vector<internal_node> nodes_;    // Indexed by node’s uid, value is node’s internal_node
  std::vector<node_id_type> i2u_;  // Indexed by node’s index, value is node’s uid

  std::vector < std::set< internal_edge > > bridge_;
};

#endif
