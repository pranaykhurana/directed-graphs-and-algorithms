#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
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


//Forward declarations for classes below so they can be used below without worrying too much about the ordering.
template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:
	
	std::vector<std::vector<bool> > adj_matrix;
	std::vector<vertex> vertices;
	
	
  //You will need to add some data members here
  //to actually represent the graph internally,
  //and keep track of whatever you need to.


public:


  directed_graph(); //A constructor for directed_graph. The graph should start empty.
  ~directed_graph(); //A destructor. Depending on how you do things, this may
  //not be necessary.
  
  bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.
  
  bool adjacent(const vertex&, const vertex&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

  void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
  void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

  void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
  void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.
  
  int index_of(const vertex&)const;
  void resize_matrix(int i);
  
  std::size_t in_degree(const vertex&) const; //Returns number of edges coming in to a vertex.
  std::size_t out_degree(const vertex&) const; //Returns the number of edges leaving a vertex.
  std::size_t degree(const vertex&) const; //Returns the degree of the vertex (both in and out edges).
  
  std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
  std::size_t num_edges() const; //Returns the total number of edges in the graph.

  std::vector<vertex> get_vertices() const; //Returns a vector containing all the vertices.
  std::vector<vertex> get_neighbours(const vertex&) const; //Returns a vector containing the neighbours of the given vertex.

  vertex_iterator<vertex> begin(); //Returns a graph_iterator pointing to the start of the vertex set.
  vertex_iterator<vertex> end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

  neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
  neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

  std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
  std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

  directed_graph<vertex> out_tree(const vertex&)const; //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
  directed_graph<vertex> in_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

  bool reachable(const vertex&, const vertex&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.

};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator {

private:
	directed_graph<vertex> current;
	size_t pos_in;


public:
  vertex_iterator(const vertex_iterator<vertex>&);
  vertex_iterator(const directed_graph<vertex>&, std::size_t);
  ~vertex_iterator();
  vertex_iterator<vertex> operator=(const vertex_iterator<vertex>&);
  bool operator==(const vertex_iterator<vertex>&) const;
  bool operator!=(const vertex_iterator<vertex>&) const;
  vertex_iterator<vertex> operator++();
  vertex_iterator<vertex> operator++(int);
  vertex operator*();
  vertex* operator->();
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator {

private:

  directed_graph<vertex> current_n;
  size_t position_n;
  vertex vert_n;

public:
  neighbour_iterator(const neighbour_iterator<vertex>&);
  neighbour_iterator(const directed_graph<vertex>&, const vertex&, std::size_t);
  ~neighbour_iterator();
  neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>&);
  bool operator==(const neighbour_iterator<vertex>&) const;
  bool operator!=(const neighbour_iterator<vertex>&) const;
  neighbour_iterator<vertex> operator++();
  neighbour_iterator<vertex> operator++(int);			
  vertex operator*();
  vertex* operator->();
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename vertex> directed_graph<vertex>::directed_graph() {}
template <typename vertex> directed_graph<vertex>::~directed_graph() {}

//this function takes in a vertex and returns the position in the list of vertices
template <typename vertex> int directed_graph<vertex>::index_of(const vertex& u) const{
	for( int i=0; i<vertices.size(); i++){
		if(vertices[i]==u)
			return i;
	}
	return -1;
}

// takes in an integer as a parameter and resizes the adj_matrix to that size
template <typename vertex> void directed_graph<vertex>::resize_matrix(int i) {
	adj_matrix.resize(i);
	for(int j=0; j<i; j++)
	{
		adj_matrix[j].resize(i);
	}
}

//searches for a given vertex in the list of vertices and returns true if found, returns false if not
template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const {
	for( int i =0; i<vertices.size(); i++){
		if(vertices[i]==u){
			return true;
		}
	}
	return false;
} 

template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) const{
	int index_u= index_of(u);
	int index_v= index_of(v);
	
	if(index_u >=0 && index_v>=0 && index_u!=index_v && index_u<vertices.size() && index_v<vertices.size())
	   {
		 return (adj_matrix[index_u][index_v]);// || adj_matrix[index_v][index_u]);
	   }
	return false;
}

// this function add a new vertex to the list of vertices and resizes the matrix after adding
template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
	vertices.push_back(u);
	resize_matrix(vertices.size());
}

//this function creates an edge between the 2 provided vertices
template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
	int index_u= index_of(u);
	int index_v= index_of(v);
	if(index_u >=0 && index_v>=0 && index_u!=index_v && index_u<vertices.size() && index_v<vertices.size() && u!=v)
		adj_matrix[index_u][index_v]=true;
	
}

// this function removes a vertex from the list of vertices and resizes the matrix after removing
template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
	if(contains(u)){
		vertices.erase(vertices.begin()+index_of(u));
	}
	resize_matrix(vertices.size());
}

// this funtion removes an edge between the vertices
template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	int index_u= index_of(u);
	int index_v= index_of(v);
	if( index_u >= 0 && index_v >= 0 && index_u != index_v && index_u < vertices.size() && index_v < vertices.size() )
		adj_matrix[index_u][index_v] = false;
}

// this functions find out the in_degree of a given vertex i.e number of edges ending on that vertex
template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) const {
	int index_u=index_of(u);
	int count = 0;
	
	for(int i = 0; i < vertices.size(); i++){
		if( adj_matrix[i][index_u] == true)
			count++;
	}
	
	return count;
}

// this function finds out the out degree of a given vertex i.e the number of edges starting on that vertex
template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u) const { 
	int index_u = index_of(u);
	int count = 0;
	
	for(int j = 0; j < vertices.size(); j++){
		if( adj_matrix[index_u][j] == true)
			count++;
	}
	
	return count;
}

//this function return the total degree of a vertex which is sum of in degree and out degree
template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u) const {
	
	return out_degree(u) + in_degree(u);
}

//this function returns the size of the list of vertices
template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const { 
	
	return vertices.size(); 
}

//this function returns the total number of edges in the graph
template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const { 
	int count = 0;
	for( int i = 0; i < vertices.size(); i++)
	{
		for( int j = 0; j < vertices.size(); j++)
		{
			if( adj_matrix[i][j] == true)
				count++;
		}
	}
	return count;
}

//this function returns the list of vertices
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() const { 
	return vertices; 
}

//this function returns the list of neighbours of a given vertex
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) const {
	int index_u = index_of(u);
	std::vector<vertex> neighbours;
	
	//out neighbours
	for(int j = 0; j < vertices.size(); j++){
		if( adj_matrix[index_u][j] == true )
			neighbours.push_back(vertices[j]);
	}
	
	return neighbours;
}

// function to set parameters for the vertex iterator to begin
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() { 
	return vertex_iterator<vertex>(*this, 0); 
}

// function to set parameters for the vertex iterator to end
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() { 
	return vertex_iterator<vertex>(*this, get_vertices().size()); 
}

// function to set parameters for the neighbour iterator to begin
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) { 
	return neighbour_iterator<vertex>(*this, u, 0); 
}

// function to set parameters for the neighbour iterator to end
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) { 
	return neighbour_iterator<vertex>(*this, u, get_neighbours(u).size()); 
}

// performs a depth first search from a given vertex and return the ordered list
template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) { 
	
	bool visited[ vertices.size() ];
	
	for (unsigned i = 0; i < vertices.size(); i++){
		visited[i] = false;
	}
	
	std::stack<vertex> unprocessed;
	
	unprocessed.push(u);
	
	std::vector<vertex> ordered;
	
	while (!unprocessed.empty()){
		
		vertex n = unprocessed.top();
		unprocessed.pop();
		
		if (!visited[index_of(n)]){
			
			visited[index_of(n)] = true;
			ordered.push_back(n);
			
			for (unsigned i = vertices.size(); i != 0; i--){
				
				if (adj_matrix[index_of(n)][i-1]){
					unprocessed.push(vertices[i-1]);
				}
			}
		}
	}
		
	return ordered;

}

// performs a breadth first search from a given vertex and return the ordered list
template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) {  
	bool visited[vertices.size()];
	
	for (unsigned i = 0; i < vertices.size(); i++){
		visited[i] = false;
	}
	
	std::queue<vertex> unprocessed;
	unprocessed.push(u);
	
	std::vector<vertex> ordered;
	
	while (!unprocessed.empty()){
		vertex n = unprocessed.front();
		unprocessed.pop();
		if (!visited[index_of(n)]){
			visited[index_of(n)] = true;
			ordered.push_back(n);
			for (unsigned i = 0; i < vertices.size(); i++){
				if (adj_matrix[index_of(n)][i]){
					unprocessed.push(vertices[i]);
				}
			}
		}
	}
		
	
	return ordered;
}

// performs created an out-edges tree from the given vertex and returns it in the form of a directed graph
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) const{
	
	int size = vertices.size();
	
	directed_graph<vertex> tree_out;
	
	bool visited[size];
	
	for (unsigned i = 0; i < size; i++){
		visited[i] = false;
	}
	
	std::queue<std::pair<vertex,vertex> > unprocessed;
	unprocessed.push({u, u});
	
	while ( !unprocessed.empty() ){
		std::pair<vertex, vertex> n = unprocessed.front();
		unprocessed.pop();
		if (!visited[index_of(n.first)]){
			visited[index_of(n.first)] = true;
			
			tree_out.add_vertex(n.first);
			tree_out.add_edge(n.second, n.first); //This only works because
											  //we check for an attempt to add
											  //a self-edge in add_edge
			for (unsigned j = 0; j < size; j++){
				if (adj_matrix[index_of(n.first)][j]){
					unprocessed.push({vertices[j], n.first});
				}
			}
			//unprocessed.pop();
		}
	}
	
	return tree_out;

}

// performs created an in-edged tree from the given vertex and returns it in the form of a directed graph
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) {

	int size = vertices.size();
	
	directed_graph<vertex> tree_in;
	
	bool visited[size];
	
	for (unsigned i = 0; i < size; i++){
		visited[i] = false;
	}
	
	std::queue<std::pair<vertex,vertex> > unprocessed;
	unprocessed.push({u, u});
	
	while ( !unprocessed.empty() ){
		std::pair<vertex, vertex> n = unprocessed.front();
		unprocessed.pop();
		if (!visited[index_of(n.first)]){
			visited[index_of(n.first)] = true;
			
			tree_in.add_vertex(n.first);
			tree_in.add_edge(n.first,n.second); //This only works because
											  //we check for an attempt to add
											  //a self-edge in add_edge
			for (unsigned j = 0; j < size; j++){
				if (adj_matrix[j][index_of(n.first)]){
					unprocessed.push({vertices[j], n.first});
				}
			}
		}
	}
	
	return tree_in;
}

// checks if a vertex v can be reached from a vertex u through a series of connected edges
template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) const {
	
	std::vector<vertex> temp;
	
	for( auto& e : out_tree(u).breadth_first(u) ){
		temp.push_back(e);
	}
	
	for( int i = 0; i < temp.size(); i++){
		if(temp[i] == v){
			return true;
		}
	}
	
	return false;
}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const vertex_iterator<vertex>& other) : current(other.current), pos_in(other.pos_in)  {}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const directed_graph<vertex>& graph, std::size_t position) : current(graph), pos_in(position) {}

template <typename vertex> vertex_iterator<vertex>::~vertex_iterator() {}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator=(const vertex_iterator<vertex>& other) { 
	current = other.current;
	pos_in = other.pos_in;
	return current, pos_in;
}

template <typename vertex> bool vertex_iterator<vertex>::operator==(const vertex_iterator<vertex>& other) const { 
	return (current.get_vertices() == other.current.get_vertices() && pos_in == other.pos_in);
}

template <typename vertex> bool vertex_iterator<vertex>::operator!=(const vertex_iterator<vertex>& other) const { 
	return ( current.get_vertices() != other.current.get_vertices() || pos_in != other.pos_in );
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++() { 
	++pos_in;
	return *this;
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++(int) { 
	auto t = this;
	++pos_in;
	return t;
}

template <typename vertex> vertex vertex_iterator<vertex>::operator*() { 
	return current.get_vertices()[pos_in]; 
}
template <typename vertex> vertex* vertex_iterator<vertex>::operator->() { 
	return current.get_vertices() + pos_in; 
}

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const neighbour_iterator<vertex>& other): current_n(other.current_n), vert_n(other.vert_n),position_n(other.position_n) {}

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position): current_n(graph), vert_n(u), position_n(position) {}

template <typename vertex> neighbour_iterator<vertex>::~neighbour_iterator() {}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator=(const neighbour_iterator<vertex>& other) { 
	
	current_n = other.current_n;
	vert_n= other.vert_n;
	position_n = other.position_n;
	
	return *this;
}

template <typename vertex> bool neighbour_iterator<vertex>::operator==(const neighbour_iterator<vertex>& other) const { 
	return (current_n.get_neighbours(vert_n) == other.current_n.get_neighbours(vert_n) && position_n == other.position_n); 
}

template <typename vertex> bool neighbour_iterator<vertex>::operator!=(const neighbour_iterator<vertex>& other) const { 
	return (current_n.get_neighbours(vert_n) != other.current_n.get_neighbours(vert_n) || position_n != other.position_n);
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++() { 
	++position_n;
	return *this;
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++(int) { 
	auto t = this;
	++position_n;
	return t;
}		

template <typename vertex> vertex neighbour_iterator<vertex>::operator*() { 
	return current_n.get_neighbours(vert_n)[position_n];
}

template <typename vertex> vertex* neighbour_iterator<vertex>::operator->() { 
	return &current_n.get_neighbours(vert_n)[position_n]; 
}


#endif