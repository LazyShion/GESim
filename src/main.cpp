/*
 * GESim
 * Last update: April 12th, 2024
 * Author: Hiroaki Shiokawa
 */
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <cassert>
#include <functional>

#include "gdb.hpp"
#include "graph_entropy.hpp"
#include "util.hpp"

using namespace std;

char *bin_file = NULL;
bool result = false;
//char *logfile = NULL;
unsigned int max_rad = 4;
unsigned int uncertainty = 1;
unsigned int gid = 0;

void
usage(char *prog_name, const string more, bool help){
  if(!help){
    cerr << "[Error] " << more << endl;
  }else{
    cerr << "[Info] " << more << endl;    
  }

  cerr << "[Usage] " << prog_name << " -i <bin_file_name> -q <query_node_id> [-r <max_radius>] [-u <uncertainty>] [-D] [-H]" << endl;
  // Display help
  if(help){
    cerr << "\t-i <sim_file_name>: Set the input file name." << endl;
    cerr << "\t-q <query_node_id>: Set a query node id." << endl;
    cerr << "\t-r <max_radius>: [Optional] Set the maximum radius." << endl;
    cerr << "\t-u <uncertainty>: [Optional] Set a uncertainty value between 0 to <max_radius>." << endl;
    cerr << "\t-D: [Optional] Display clusteering results" << endl;
    cerr << "\t-H: [Optional] Display help menu." << endl;
    cerr << "" << endl;
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
  if(argc <= 1){
    usage(argv[0], "Arguments are missing.", true);    
  }
  
  for(int i = 1; i < argc; i++){
    if(argv[i][0] == '-'){
      switch(argv[i][1]){
      case 'i':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	bin_file = argv[i+1];
	i+=1;
	break;
      case 'q':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	gid = (unsigned int) atoi(argv[i+1]);
	i+=1;
	break;
      case 'r':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	max_rad = (unsigned int) atoi(argv[i+1]);
	i+=1;
	break;
      case 'u':
	check_args(i+1, argc, argv[0], "Invalid arguments at "+string(argv[i]));
	uncertainty = (unsigned int) atoi(argv[i+1]);
	i+=1;
	break;	
      case 'D':
	result = true;
	break;
      case 'H':
	usage(argv[0], "Help menu", true);
	break;
      default:
	usage(argv[0], "Unknown option", false);
      }
    }
  }

  if(bin_file == NULL){
    usage(argv[0], "<bin_file_name> is missing.", false);
  }

}


int
main(int argc, char **argv){
  parse_args(argc, argv);
  GraphDB gdb(bin_file);
  GraphEntropy g_entropy(&gdb, max_rad, uncertainty);  

  //vector<double> results = g_entropy.search_all(gid);

  if(result){
    // For All pairs sim test on zinc1000
    /*
    for(unsigned int gid = 0, end=gdb.N; gid<end; ++gid){
      for(unsigned int gid1 = 0, end1=gdb.N; gid1<end1; ++gid1){
	double sim = g_entropy.comp_QJS(gid, gid1);
	sim = 1 - sim;
	//cout << gid << "\t" <<gid1 << "\t" << fixed << setprecision(5) << sim << endl;
	cout << fixed << setprecision(5) << sim << endl;
      }
    }
    */
    // For default query test
    vector<double> result_all = g_entropy.graph_entropy_all(615);
    vector<double> result_all_best = g_entropy.graph_entropy_all_best(615);

    for(unsigned int i=0; i<result_all.size(); ++i){
      cout << "615\t" << i << "\t" << result_all[i] << "\t" << result_all_best[i] << endl;
    }
    
    /*
    for(unsigned int gid1 = 0, end1=gdb.N; gid1<end1; ++gid1){
      double sim = g_entropy.comp_QJS_Best(gid, gid1);
      //double sim = g_entropy.comp_QJS(gid, gid1);
      cout << gid << "\t" <<gid1 << "\t" << fixed << setprecision(5) << sim << endl;
    }
    */
    /*
    cout << "Original: " << 1 - g_entropy.comp_QJS(0, 1) << endl;
    cout << "Updated: " << 1 - g_entropy.comp_QJS_Best(0, 1) << endl;
    */
  }
  
}

