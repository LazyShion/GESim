//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include "graph_entropy.hpp"

using namespace std;

//namespace py = pybind11;

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

// Test codes
vector<pair<double, unsigned int>>
GraphEntropy::search_all(unsigned int query){
  vector<pair<double, unsigned int>> result(gdb->N);
  
  for(unsigned int gid=0, N=(*gdb).N; gid<N; ++gid){
    result[gid] = make_pair(comp_QJS(query, gid), gid);
  }

  return result;
}

// Return a graph entropy between graphs gid1 and gid2
double
GraphEntropy::graph_entropy(unsigned int gid1, unsigned int gid2){
  return comp_QJS(gid1, gid2);
}

// Return match mapping between g1 and g2
vector<vector<int>>
GraphEntropy::match_mapping(unsigned int gid1, unsigned int gid2){
  pair<vector<int>, vector<int>> p = align_match(gid1, gid2);
  vector<vector<int>> mapping(2);
  mapping[0] = p.first;
  mapping[1] = p.second;

  return mapping;
}

// Return graph entropies between graph gid1 and each of gdb.
vector<double>
GraphEntropy::graph_entropy_all(unsigned int gid1){
  vector<double> result(gdb->N-1);
  unsigned int i=0;

  for(unsigned int gid2=0, N=(*gdb).N; gid2<N; ++gid2){
    if(gid2 != gid1){
      result[i] = comp_QJS(gid1, gid2);
      ++i;
    }
  }
  
  return result;
}
/*
PYBIND11_MODULE(graph_entropy, m) {
    py::class_<GraphEntropy>(m, "GraphEntropy")
        .def(py::init<>())
        .def(py::init<GraphDB*, unsigned int>())
        .def("set_query", &GraphEntropy::set_query)
        .def("search_all", &GraphEntropy::search_all)
        .def("graph_entropy", &GraphEntropy::graph_entropy)
        .def("graph_entropy_all", &GraphEntropy::graph_entropy_all)
        .def("match_mapping", &GraphEntropy::match_mapping);
}
*/
