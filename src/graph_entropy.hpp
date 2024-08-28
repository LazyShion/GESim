#ifndef GraphEntropy_H
#define GraphEntropy_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <functional>
#include <bitset>
#include <utility>
#include <queue>

#include "util.hpp"
#include "gdb.hpp"

using namespace std;

class GraphEntropy{
private:
  // ---------------------
  // Members
  // ---------------------
  GraphDB *gdb; // graph database
  vector<unsigned int> *q_nodes;   // query_nodes
  vector<unsigned int> *q_labels;  // query_labels
  vector<unsigned int> *q_edges;   // query_edges
  vector<unsigned int> *q_weights; // query_weights
  const static unsigned int FP_LEN = 1024; // length of a finger print
  vector<bitset<FP_LEN>> ecfp; // finger print vector
  unsigned int r; // max_rad
  
public:
  // ---------------------
  // Members
  // ---------------------
  
  
  // ---------------------
  // Constructor
  // ---------------------
  
  GraphEntropy();
  GraphEntropy(GraphDB *db, unsigned int max_rad);

  
  // ---------------------
  // Functions
  // ---------------------
  void set_query(vector<unsigned int> *nodes, vector<unsigned int> *labels, vector<unsigned int> *edges, vector<unsigned int> *weights);
  vector<pair<double, unsigned int>> search_all(unsigned int query);
  double graph_entropy(unsigned int gid1, unsigned int gid2);
  vector<double> graph_entropy_all(unsigned int gid1);
  
  // ---------------------
  // Inlines 
  // ---------------------

  // compute struct infomation for graph gid
  inline double
  comp_SI(unsigned int gid){
    double si = 0.0;
    unsigned int volG = (*gdb).num_edges_gid(gid);
    unsigned int num_nodes = (*gdb).num_nodes_gid(gid);

    for(unsigned int i=0; i<num_nodes; ++i){
      unsigned int deg = (*gdb).degree_gid(gid, i);
      if(deg>0){
	double x = (double)deg/volG;
	si += x*log(x);
      }
    }
    return -si;
  }

  // Graph alignment between g1 and g2 by GRAAL
  inline double
  align_graphs(unsigned int g1, unsigned int g2){
    // Initialization
    unsigned int num_nodes_g1 = gdb->num_nodes_gid(g1);
    unsigned int num_nodes_g2 = gdb->num_nodes_gid(g2);
    if(num_nodes_g1 > num_nodes_g2){
      swap(g1, g2);
      swap(num_nodes_g1, num_nodes_g2);
    }

    // Costs computation
    vector<double> costs(num_nodes_g1*num_nodes_g2, -1);
    vector<unsigned int>::iterator node_g1 = gdb->nodes_gid(g1);
    vector<unsigned int>::iterator node_g2 = gdb->nodes_gid(g2);
    unsigned int prefix_g1 = node_g1 - gdb->nodes.begin();
    unsigned int prefix_g2 = node_g2 - gdb->nodes.begin();
    unsigned int vol_g1 = gdb->num_edges_gid(g1);
    unsigned int vol_g2 = gdb->num_edges_gid(g2);
    
    for(unsigned int i=0; i<num_nodes_g1; ++i){
      for(unsigned int j=0; j<num_nodes_g2; ++j){
	bitset<FP_LEN> bs1 = ecfp[prefix_g1+i] & ecfp[prefix_g2+j];
	double Tsim = (double)bs1.count();
	costs[num_nodes_g2*i+j] = Tsim;
      }
    }

    // Alignment by GRAAL
    vector<int> align_g1(num_nodes_g1, -1);
    vector<int> align_g2(num_nodes_g2, -1);
    GRAAL(g1, g2, align_g1.begin(), align_g1.size(), align_g2.begin(), align_g2.size(), costs.begin());    
    
    // Merge g1 and g2 based on the alignment
    vector<unsigned int> degrees_g1(num_nodes_g1, 0);
    vector<unsigned int> degrees_g2(num_nodes_g2, 0);
    unsigned int total_deg = 0;

    for(unsigned int i=0; i<num_nodes_g1; ++i){
      int pos = align_g1[i];

      if(pos != -1 && costs[num_nodes_g2*i + pos] == 1){
	degrees_g1[i] = gdb->merged_degree(g1, i, g2, pos);
      }else{
	degrees_g1[i] = gdb->degree_gid(g1, i);
      }
      total_deg += degrees_g1[i];
    }
    
    for(unsigned int j=0; j<num_nodes_g2; ++j){
      int pos = align_g2[j];
      if(pos == -1){
	degrees_g2[j] = gdb->degree_gid(g2, j);
	total_deg += degrees_g2[j];
      }      
    }
    
    // Compute graph entropy (structural information)
    double gent = 0.0;
    for(unsigned int i=0; i<num_nodes_g1; ++i){
      unsigned int deg = degrees_g1[i];
      int pos = align_g1[i];
      double cost = r;

      if(pos != -1){
	cost = cost - costs[num_nodes_g2*i+pos];
      }
      
      if(deg > 0){	
	double x = (double)deg/total_deg;
	gent += x * log(x);
      }
    }
    
    for(unsigned int j=0; j<num_nodes_g2; ++j){
      unsigned int deg = degrees_g2[j];
      int pos = align_g2[j];
      double cost = r;

      if(pos != -1){
	cost = cost - costs[num_nodes_g2*pos+j];
      }
      
      if(deg > 0){
	double x = (double)deg/total_deg;
	gent += x * log(x);
      }
    }

    return -gent;
  }


  // compute QJS distance between graphs g1 and g2
  inline double
  comp_QJS(unsigned int g1, unsigned int g2){
    double diff = align_graphs(g1, g2) - (comp_SI(g1)+comp_SI(g2))/2;

    if(diff<=0){
      return 0;
    }else{
      return sqrt(diff);
    }
  }

  // generate ecfp code for a node nid in a graph gid
  inline bitset<FP_LEN>
  generate_ecfp_nid(GraphDB *gdb, unsigned int gid, unsigned int nid, unsigned int max_rad){
    bitset<FP_LEN> fp(0);
    
    for(unsigned int dir=0; dir<=max_rad; ++dir){
      string code = generate_PATH_FP(gdb, gid, nid, dir);
      unsigned int pos = hash<std::string>()(code)%FP_LEN;

      // if pos is conflicted, ++pos until fp[pos] is satisfied.
      while(true){
	if(!fp[pos]){
	  fp.set(pos);
	  break;
	}else{
	  pos++;
	  pos = pos%FP_LEN;
	}
      }      
    }
    
    return fp;
  }

  // generate ecfp code for a graph gid
  inline bitset<FP_LEN>
  generate_ecfp_gid(GraphDB *gdb, unsigned int gid, unsigned int max_rad){
    bitset<FP_LEN> fp(0);
    
    for(unsigned int nid=0, end_nid=gdb->num_nodes_gid(gid); nid<end_nid; ++nid){
      fp |= generate_ecfp_nid(gdb, gid, nid, max_rad);
    }
    
    return fp;
  }
  
  // Generate Path-based FP for node nid in graph gid with maximum radius (max_rad)
  inline string
  generate_PATH_FP(GraphDB *gdb, unsigned int gid, unsigned int nid, unsigned int max_rad){    
    list<string> PATHS;
    list<pair<unsigned int, pair<unsigned int, pair<int, string>>>> queue;
    queue.push_back(make_pair(max_rad, make_pair(nid, make_pair(-1, "#"))));

    while(!queue.empty()){
      pair<unsigned int, pair<unsigned int, pair<int, string>>> elm  = queue.front();
      queue.pop_front();
      unsigned int rad = elm.first;
      unsigned int node = (elm.second).first;
      unsigned int label = (*gdb).node_label(gid, node);
      int parent = ((elm.second).second).first;
      vector<unsigned int>::iterator it = (*gdb).neighbors(gid, node);
      vector<double>::iterator it_w = (*gdb).neighbor_weights(gid, node);
      unsigned int degree = (*gdb).degree_gid(gid, node);
      string path = ((elm.second).second).second+"["+to_string(label)+":"+to_string(degree)+"]";

      if(rad !=max_rad && degree == 1){
	PATHS.push_back(path);
      }else{
	for(unsigned int i=0; i<degree; ++i){
	  unsigned int neigh = *(it+i);
	  unsigned int weight = *(it_w+i)*2;
	  string sub_path = "-("+to_string(weight)+")-";
	  
	  if(neigh != parent){
	    if(rad>1){
	      queue.push_back(make_pair(rad-1, make_pair(neigh, make_pair((int)node, path+sub_path))));
	    }else{
	      PATHS.push_back(path+sub_path);
	    }
	  }
	}
      }
    }

    string code = "";
    PATHS.sort();
    list<string>::iterator it=PATHS.begin(), end=PATHS.end();
    
    while(it != end){
      code += *it;
      ++it;
    }
    
    return code;
  }

  // Run GRAAL old between g1 and g2
  inline void
  GRAAL(unsigned int g1, unsigned int g2,
	vector<int>::iterator it_g1, unsigned int size_g1,
	vector<int>::iterator it_g2, unsigned int size_g2,
	vector<double>::iterator it_costs){

    for(unsigned int val=r; val>0; --val){
      
      for(unsigned int i=0; i<size_g1; ++i){
	if(*(it_g1+i)!=-1){
	  continue;
	}
	
	for(unsigned int j=0; j<size_g2; ++j){
	  unsigned int index = size_g2*i+j;
	  if(*(it_g2+j) == -1){
	    if(*(it_costs+index) > r-1){
	      *(it_g1 + i) = j;
	      *(it_g2 + j) = i;
	      break;
	    }
	  }
	}
      }      
    }    
  }  
  
  
};

#endif //GraphEntropy_H
