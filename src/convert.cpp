#include <pybind11/pybind11.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <cassert>

#include "util.hpp"

using namespace std;

namespace py = pybind11;

char *input_file = NULL;
char *output_file = NULL;
unsigned int N=0; // # of graphs
unsigned int V=0; // # of total nodes
unsigned int E=0; // # of total edges
unsigned int L=0; // # of label types

struct Graph{
  unsigned int gid;
  vector<string> nodes; //nid:nlabel
  vector<vector<pair<unsigned int, double>>> edges; // src: vec(<dst, weight>)
};

void
usage(char *prog_name, const string more, bool help){
  cerr << "[Error] " << more << endl;
  cerr << "[Usage] " << prog_name << " -i <input_file_name> -o <output_file_name> [-h]" << endl;
  // Display help
  if(help){
    cerr << "\t-h: Display help menu" << endl;
    cerr << "\t-i <input_file_name>: Set the input file name." << endl;
    cerr << "\t-o <output_file_name>: Set the output file name." << endl;
  }
  cerr << endl;
  exit(0);
}

void
check_args(int pos, int argc, char *prog_name, string more){
  if(pos >= argc){
    usage(prog_name, more, false);
  }
}

void
parse_args(int argc, char **argv) {
  for(int i = 1; i < argc; i++){
    if(argv[i][0] == '-'){
      switch(argv[i][1]){
      case 'i':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	input_file = argv[i+1];
	i++;
	break;
      case 'o':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	output_file = argv[i+1];
	i++;
	break;
      case 'h':
	usage(argv[0], "Help menu", true);
	break;
      default:
	usage(argv[0], "Unknown option", false);
      }
    }
  }

  if(input_file == NULL){
    usage(argv[0], "<input_file_name> is missing.", false);
  }
  
  if(output_file == NULL){
    usage(argv[0], "<output_file_name> is missing.", false);
  }
}

void
read_file_impl(ifstream &finput,
    vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){

  string buf;
  unsigned int gid;
  unsigned int nid;
  string nlabel;
  unsigned int eid, src, dst;
  double weight;

  vector<Graph> gdb(0);
  Graph g = {};
  N=0; // # of graphs
  V=0; // # of total nodes
  E=0; // # of total edges
  L=0; // # of label types  
  unordered_map<string, unsigned int> l2i; // mapping a label to an integer

  
  while(finput >> buf){

    if(buf == "t"){
      finput >> gid;
      g.gid = gid;
      gdb.push_back(g);
      ++N;
    }else if(buf == "v"){
      finput >> nid >> nlabel;
      if(gdb[N-1].nodes.size()<=nid){
	gdb[N-1].nodes.resize(nid+1);
      }
      
      gdb[N-1].nodes[nid]= nlabel;
      ++V;
      
      if(l2i[nlabel] == 0){
        l2i[nlabel] = (unsigned int)stoi(nlabel);
        //l2i[nlabel] = L;
	//      ++L;
      }
      
    }else if(buf == "e"){
      finput >> eid >> src >> dst >> weight;
      if(gdb[N-1].edges.size()==0){
	gdb[N-1].edges.resize(gdb[N-1].nodes.size());
      }
      gdb[N-1].edges[src].push_back(make_pair(dst, weight));
      gdb[N-1].edges[dst].push_back(make_pair(src, weight));
      E+=2;
    }
  }
  
  // construct index
  (*graph_index).resize(N, 0);
  (*node_index).resize(V, 0);
  (*label_index).resize(V, 0);  
  (*edge_index).resize(E, 0);
  (*weight_index).resize(E, 0);
  unsigned int prev=0;

  for(unsigned int i=0; i<N; ++i){
    // graph index construction
    if(i > 0){
      prev = (*graph_index)[i-1];
    }
    (*graph_index)[i] = prev + gdb[i].nodes.size();

    for(unsigned int j=0, end_j = gdb[i].nodes.size(); j<end_j; ++j){
      // node index construction
      if(prev+j == 0){
	(*node_index)[0] = gdb[i].edges[j].size();
	(*label_index)[0] = l2i[gdb[i].nodes[j]];
      }else{
	(*node_index)[prev+j] = (*node_index)[prev+j-1] + gdb[i].edges[j].size();
	(*label_index)[prev+j] = l2i[gdb[i].nodes[j]];
      }

      // edge index construction
      for(unsigned int k=0, end_k=gdb[i].edges[j].size(); k<end_k; ++k){
	if(prev+j == 0){
	  (*edge_index)[k] = gdb[i].edges[j][k].first;
	  (*weight_index)[k] = gdb[i].edges[j][k].second;
	}else{
	  (*edge_index)[(*node_index)[prev+j-1]+k] = gdb[i].edges[j][k].first;
	  (*weight_index)[(*node_index)[prev+j-1]+k] = gdb[i].edges[j][k].second;
	}
      }
    }    
  }
  
}

void
read_file(vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){

  // read input file
  ifstream finput;
  finput.open(input_file, fstream::in);
  read_file_impl(finput, graph_index, node_index, label_index, edge_index, weight_index);
  finput.close();
}

void
read_file(const string &input_file_name,
    vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){

  // read input file
  ifstream finput;
  finput.open(input_file_name, fstream::in);
  read_file_impl(finput, graph_index, node_index, label_index, edge_index, weight_index);
  finput.close();
}

void
write_binary_impl(ofstream &foutput,
    vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){
  
  // output statistics
  foutput.write((char *)(&N), 4);
  foutput.write((char *)(&V), 4);
  foutput.write((char *)(&E), 4);

  // output graph_index
  for(unsigned int i=0; i<N; ++i){
    unsigned int dump = (*graph_index)[i];
    foutput.write((char *)(&dump), 4);
  }
  
  // output node_index
  for(unsigned int i=0; i<V; ++i){
    unsigned int dump = (*node_index)[i];
    foutput.write((char *)(&dump), 4);
  }

  // output label_index
  for(unsigned int i=0; i<V; ++i){
    unsigned int dump = (*label_index)[i];
    foutput.write((char *)(&dump), 4);
  }

  // output edge_index
  for(unsigned int i=0; i<E; ++i){
    unsigned int dump = (*edge_index)[i];
    foutput.write((char *)(&dump), 4);
  }

  // output weight_index
  for(unsigned int i=0; i<E; ++i){
    double dump = (*weight_index)[i];
    foutput.write((char *)(&dump), 8);
  }
}

void
write_binary(vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){
  
  ofstream foutput;
  foutput.open(output_file, fstream::out | fstream::binary);
  write_binary_impl(foutput, graph_index, node_index, label_index, edge_index, weight_index);
  foutput.close();
}

void
write_binary(const string &output_file_name,
    vector<unsigned int> *graph_index,
	  vector<unsigned int> *node_index,
	  vector<unsigned int> *label_index,
	  vector<unsigned int> *edge_index,
	  vector<double> *weight_index){
  
  ofstream foutput;
  foutput.open(output_file_name, fstream::out | fstream::binary);
  write_binary_impl(foutput, graph_index, node_index, label_index, edge_index, weight_index);
  foutput.close();
}

void
convert_graph_to_binary(const string &input_file_name, const string &output_file_name){
  vector<unsigned int> graph_index(0);
  vector<unsigned int> node_index(0);
  vector<unsigned int> label_index(0);
  vector<unsigned int> edge_index(0);
  vector<double> weight_index(0);

  read_file(input_file_name, &graph_index, &node_index, &label_index, &edge_index, &weight_index);
  write_binary(output_file_name, &graph_index, &node_index, &label_index, &edge_index, &weight_index);
}

PYBIND11_MODULE(convert, m) {
  m.def("convert_graph_to_binary", &convert_graph_to_binary, "Convert graph data from text to binary format",
      py::arg("input_file_name"), py::arg("output_file_name"));
}

int
main(int argc, char **argv){
  if(argc > 1){
    parse_args(argc, argv);
  }
  // -----------------------------------------------------------------------------------
  // Main procedures
  // -----------------------------------------------------------------------------------

  vector<unsigned int> graph_index(0);
  vector<unsigned int> node_index(0);
  vector<unsigned int> label_index(0);
  vector<unsigned int> edge_index(0);
  vector<double> weight_index(0);
  
  // Read input file
  read_file(&graph_index, &node_index, &label_index, &edge_index, &weight_index);
  
  // Write output file
  write_binary(&graph_index, &node_index, &label_index, &edge_index, &weight_index);

  
  
}
