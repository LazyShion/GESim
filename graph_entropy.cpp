#include "graph_entropy.hpp"

using namespace std;

GraphEntropy::GraphEntropy(){}


GraphEntropy::GraphEntropy(GraphDB *db, unsigned int max_rad){
  gdb = db;
  ecfp.resize(gdb->V);
  unsigned int pos = 0;
  r = max_rad;

  for(unsigned int gid=0; gid<gdb->N; ++gid){
    for(unsigned int nid=0, deg=gdb->num_nodes_gid(gid); nid<deg; ++nid){
      ecfp[pos] = generate_ecfp_nid(gdb, gid, nid, max_rad);
      ++pos;
    }
  }
  
}


void
GraphEntropy::set_query(vector<unsigned int> *nodes,
			vector<unsigned int> *labels,
			vector<unsigned int> *edges,
			vector<unsigned int> *weights){
  q_nodes = nodes;
  q_labels = labels;
  q_edges = edges;
  q_weights = weights;

}


vector<pair<double, unsigned int>>
GraphEntropy::search_all(unsigned int query){
  vector<pair<double, unsigned int>> result(gdb->N);
  
  for(unsigned int gid=0, N=(*gdb).N; gid<N; ++gid){
    result[gid] = make_pair(comp_QJS(query, gid), gid);
  }

  return result;
}


