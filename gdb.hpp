#ifndef GraphDB_H
#define GraphDB_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <functional>
#include <list>

#include "util.hpp"

using namespace std;

class GraphDB{
private:
  // ---------------------
  // Members
  // ---------------------


public:
  // ---------------------
  // Members
  // ---------------------
  unsigned int N; // # of graphs
  unsigned int V; // # of total nodes
  unsigned int E; // # of total edges
  vector<unsigned int> graphs; // graph_index
  vector<unsigned int> nodes; // node_index
  vector<unsigned int> labels; // label_index
  vector<unsigned int> edges; // edge_index
  vector<double> weights; // weight_index
  
  // ---------------------
  // Constructor
  // ---------------------

  GraphDB();
  GraphDB(char *bin_file);
  
  // ---------------------
  // Functions
  // ---------------------

  void read_from_binary(char *bin_file);

  // ---------------------
  // Inlines 
  // ---------------------

  // Return # of nodes in graph gid
  inline unsigned int
  num_nodes_gid(unsigned int gid){
    if(gid == 0){
      return graphs[gid];
    }else{
      return graphs[gid] - graphs[gid-1];
    }
  }
  
  // Return # of edges in graph gid
  inline unsigned int
  num_edges_gid(unsigned int gid){
    if(gid == 0){
      return nodes[graphs[gid]-1];
    }else{
      return nodes[graphs[gid]-1] - nodes[graphs[gid-1]-1];
    }
  }

  // Return a label of node u in graph gid
  inline unsigned int
  node_label(unsigned int gid, unsigned int u){
    if(gid == 0){
      return labels[u];
    }else{
      return labels[graphs[gid-1] + u];
    }
  }

  // Return a head iterator of nodes in graph gid
  vector<unsigned int>::iterator
  nodes_gid(unsigned int gid){
    vector<unsigned int>::iterator it = nodes.begin();
    if(gid == 0){
      return it;
    }else{
      return it + graphs[gid-1];
    }
  }

  // Return a head iterator of edges in graph gid
  vector<unsigned int>::iterator
  edges_gid(unsigned int gid){
    vector<unsigned int>::iterator it = edges.begin();
    if(gid == 0){
      return it;
    }else{
      return it + nodes[graphs[gid-1]-1];
    }
  }
  
  // Return degree of node u in graph gid
  inline unsigned int
  degree_gid(unsigned int gid, unsigned int u){
    unsigned int prefix = 0;
    if(gid == 0){
      if(u == 0){
	return nodes[u];
      }else{
	return nodes[u] - nodes[u-1];
      }
    }else{
      prefix = graphs[gid-1];
      return (nodes[prefix + u]) - (nodes[prefix + u -1]);
    }
  }

  // Return a head iterator of neighbor nodes
  inline vector<unsigned int>::iterator
  neighbors(unsigned int gid, unsigned int u){
    vector<unsigned int>::iterator it = edges.begin();
    unsigned int prefix = 0;
    
    if(gid == 0){
      if(u == 0){
	return it;
      }else{
	return it + nodes[u-1];
      }
    }else{
      prefix = graphs[gid - 1];
      return it + (nodes[prefix + u -1]);
    }
  }

  // Return a head iterator of neighbor nodes' weight
  inline vector<double>::iterator
  neighbor_weights(unsigned int gid, unsigned int u){
    vector<double>::iterator it = weights.begin();
    unsigned int prefix = 0;

    if(gid == 0){
      if(u == 0){
	return it;
      }else{
	return it + nodes[u-1];
      }
    }else{
      prefix = graphs[gid - 1];
      return it + (nodes[prefix + u - 1]);
    }
  }

  // Return a merged degree
  inline unsigned int
  merged_degree(unsigned int g1, unsigned int pos1, unsigned int g2, unsigned int pos2){
    vector<unsigned int>::iterator neigh_g1 = neighbors(g1, pos1);
    vector<unsigned int>::iterator neigh_g2 = neighbors(g2, pos2);
    vector<double>::iterator weights_g1 = neighbor_weights(g1, pos1);
    vector<double>::iterator weights_g2 = neighbor_weights(g2, pos2);
    unsigned int degree_g1 = degree_gid(g1, pos1);
    unsigned int degree_g2 = degree_gid(g2, pos2);
    unsigned int num_matches = 0;

    
    for(unsigned int i=0; i<degree_g1; ++i){
      for(unsigned int j=0; j<degree_g2; ++j){

	// match
	if((*(neigh_g1+i) == *(neigh_g2+j)) && (*(weights_g1+i) == *(weights_g2+j))){
	  ++num_matches;
	  break;
	}
      }
    }
    
    return degree_g1 + degree_g2 - num_matches;
  }

  
};

#endif //GraphDB_H
