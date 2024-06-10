#include <pybind11/pybind11.h>
#include "gdb.hpp"

using namespace std;
namespace py = pybind11;

GraphDB::GraphDB(){}

GraphDB::GraphDB(char *bin_file){
  read_from_binary(bin_file);
}

void
GraphDB::read_from_binary(char *bin_file){
  ifstream finput;
  finput.open(bin_file, fstream::in | fstream::binary);

  finput.read((char *)(&N), 4);
  finput.read((char *)(&V), 4);
  finput.read((char *)(&E), 4);

  graphs.resize(N);
  nodes.resize(V);
  labels.resize(V);
  edges.resize(E);
  weights.resize(E);

  finput.read((char *)(&graphs[0]),  N*4);
  finput.read((char *)(&nodes[0]),   V*4);
  finput.read((char *)(&labels[0]),  V*4);
  finput.read((char *)(&edges[0]),   E*4);
  finput.read((char *)(&weights[0]), E*8);
  
  finput.close();  
}

PYBIND11_MODULE(gdb, m) {
    py::class_<GraphDB>(m, "GraphDB")
        .def(py::init<>())
        .def(py::init<char*>())
        .def("read_from_binary", &GraphDB::read_from_binary)
        .def("num_nodes_gid", &GraphDB::num_nodes_gid)
        .def("num_edges_gid", &GraphDB::num_edges_gid)
        .def("node_label", &GraphDB::node_label)
        .def("nodes_gid", &GraphDB::nodes_gid)
        .def("edges_gid", &GraphDB::edges_gid)
        .def("degree_gid", &GraphDB::degree_gid)
        .def("neighbors", &GraphDB::neighbors)
        .def("neighbor_weights", &GraphDB::neighbor_weights)
        .def("merged_degree", &GraphDB::merged_degree)
        // Expose the public members
        .def_readwrite("N", &GraphDB::N)
        .def_readwrite("V", &GraphDB::V)
        .def_readwrite("E", &GraphDB::E)
        .def_readwrite("graphs", &GraphDB::graphs)
        .def_readwrite("nodes", &GraphDB::nodes)
        .def_readwrite("labels", &GraphDB::labels)
        .def_readwrite("edges", &GraphDB::edges)
        .def_readwrite("weights", &GraphDB::weights);
}

