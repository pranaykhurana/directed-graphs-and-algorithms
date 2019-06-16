/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>
#include <iostream>

#include "directed_graph.hpp"

/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */

//performs  a topological sort and returns true/false based on if there are any edges remaining in the graph
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  directed_graph<vertex> g(d);
  std::list<vertex> topo_order;
  std::set<vertex> s;
  
  for ( auto e = g.begin(); e != g.end(); e++ ){
    if( g.in_degree(*e) == 0 )
      s.insert(*e);
  }
  
  while( !s.empty() ){
    auto v = *s.begin();
    s.erase(v);
    topo_order.push_back(v);
    for ( auto a : g ) if ( g.adjacent(v, a) ) {
      g.remove_edge(v,a);
      if( g.in_degree(a) == 0 ) {
        s.insert(a);
      }
    }
  } 
  if( g.num_edges() > 0)
    return false; 
  return true; 
}


/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
 
// takes in a graph as a parameter and performs topological sorting on the graph based on kahn's algorithm and returns the ordered list
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
  directed_graph<vertex> g(d);
  std::list<vertex> topo_order;
  std::set<vertex> s;
  
  for ( auto e = g.begin(); e != g.end(); e++ ){
    if( g.in_degree(*e) == 0 )
      s.insert(*e);
  }
  
  while( !s.empty() ){
    auto v = *s.begin();
    s.erase(v);
    topo_order.push_back(v);
    for ( auto a : g ){ 
      if ( g.adjacent(v, a) ) {
        g.remove_edge(v,a);
        if( g.in_degree(a) == 0 ) {
          s.insert(a);
        }
      }
    }  
  }
  
  return topo_order;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
 
//performs a topological sort(kahn's algorithm) on the graph and makes sure every vertex is visited once and return true/false
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
  std::list<vertex> ordered = topological_sort(d);
  
  while( ordered.size() > 1 )
  {
    vertex v = ordered.front();
    ordered.erase(ordered.begin());
    if( !d.adjacent(v,ordered.front()) )
      return false; 
  }
  
  return true;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
 // takes in the graph as a parameter and performs a breadth first search, checks for the adjacency between different vertices and returns a ordered list of weakly connected components
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
   std::vector<std::vector<vertex>> comp;
   std::unordered_map<vertex, bool> visited;
   std::stack<vertex> unprocessed;
   std::vector<vertex> ordered;
   
   for( auto i = d.begin(); i != d.end(); i++){
     visited[*i] = false;
   }
   
   for( auto a = d.begin(); a != d.end(); a++ ){
     unprocessed.push(*a);
      while( !unprocessed.empty() ){
        vertex n = unprocessed.top();
        unprocessed.pop();
        if( visited[n] == false ){
          ordered.push_back(n);
          visited[n] = true;
          for( auto x = d.begin(); x != d.end(); x++ ){
            if( d.adjacent(n,*x) || d.adjacent(*x,n) )
            unprocessed.push(*x);
          }
        }
      }
     if( ordered.size() != 0 ){
     comp.push_back(ordered);
     ordered.clear();
     } 
   } 
   
  return comp;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
//part -1 // takes in the graph and sets up all the parameters and sends it to part 2
template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
   std::vector<std::vector<vertex>> strong_comp;
   int count = 0;
   std::stack<vertex> s;
   std::unordered_map<vertex, int> index;
   std::unordered_map<vertex, int> low;
   std::unordered_map<vertex, bool> onS;
   
   for( auto a = d.begin(); a != d.end(); a++ ){
     onS[*a] = false;
     index[*a] = -1;
   }
   
   for( auto a = d.begin(); a != d.end(); a++ ){
     if( index[*a] == -1 ){
       strongconnect(*a, d, count, s, index, low, onS, strong_comp);
     }
   } 
   
 return strong_comp;
}

// performs tarjans algorithm on by using the provided pointers from part 1 to return a vector matrix containing the strongly connected components
template<typename vertex>
void strongconnect(const vertex & u, const directed_graph<vertex> & d, int & count, std::stack<vertex> & s,std::unordered_map<vertex, int> & index, std::unordered_map<vertex, int> & low, std::unordered_map<vertex, bool> & onS, std::vector<std::vector<vertex>> & strong_comp){
  s.push(u);
  onS[u] = true;
  index[u] = count;
  low[u] = count;
  count = count + 1;
  
  for( auto n = d.nbegin(u); n != d.nend(u); n++ ){
      if( index[*n] == -1 ){
        vertex v = *n;
        strongconnect(v, d, count, s, index, low, onS, strong_comp);
        low[u] = std::min(low[u], low[v]);
      }
      if( onS[*n] == true ){
        low[u] =  std::min(low[u], index[*n]);
      }
  }
  
  if( low[u] == index[u] ){
    std::vector<vertex> new_component;
    vertex x;
    do {
      x = s.top();
      s.pop();
      onS[x] = false;
      new_component.push_back(x);
    }
    while( x != u );
    strong_comp.push_back(new_component);
  }
} 

/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */

// calculated the shortest path from a given vertex to every other vertex by performing a breadth first search and returns a map containing the distances and vertices
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  std::unordered_map<vertex, std::size_t> distances;
  std::size_t dist = 1;
  std::unordered_set<vertex> visited;
  
  for( auto v = d.begin(); v != d.end(); v++ ){
      distances[*v] = d.num_vertices() + 1;
  }
  distances[u] = 0;
  std::queue<vertex> unprocessed;
	unprocessed.push(u);
  while ( !unprocessed.empty() ){
		vertex n = unprocessed.front();
		unprocessed.pop();
    if( visited.count(n) == 0 ){
      for( auto i = d.nbegin(n); i != d.nend(n); i++ ){
        unprocessed.push(*i);
        if( distances[*i] > distances[n] )
          distances[*i] = distances[n] + 1;
      }
      dist++;
    }
    visited.insert(n);
  } 
  return distances;
}

int main() {

	
}