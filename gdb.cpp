#include "gdb.hpp"

using namespace std;

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
