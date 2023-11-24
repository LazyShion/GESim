/*
 * My util functions (version 0.2.0)
 * Last update: April 22, 2020
 * Author: Hiroaki Shiokawa
 */
#ifndef UTIL_H
#define UTIL_H

#include<stdlib.h>
#include<fstream>
#include<chrono>
#include<unordered_map>
#include<list>
#include<random>
#include<functional>

using namespace std;

static string LOGS;
static unordered_map<string,string> LOG_MAP;

typedef chrono::system_clock::time_point Chrono_t;

// macro for debug
#define debug(var) do{cout << #var << ": "; view(var); }while(0)
template<typename T> void view(T e){cout << e << endl;}
template<typename T> void view(const vector<T>& v){for(const auto& e:v){cout << e << " ";} cout << endl;}
template<typename T> void view(const vector<vector<T>>& vv){for(const auto& v:vv){view(v);}}

// get a time point
static Chrono_t get_time(){
  return chrono::system_clock::now();
}

// get a running time in milliseconds
static double get_time_millisec(const Chrono_t *start, const Chrono_t *end){
  return static_cast<double>(chrono::duration_cast<chrono::microseconds>(*end - *start).count()/ 1000.0);
}

// get a running time in seconds
static double get_time_sec(const Chrono_t *start, const Chrono_t *end){
  return static_cast<double>(chrono::duration_cast<chrono::microseconds>(*end - *start).count()/ 1000000.0);
}

// output logs in JSON format
static void output_json_logs(const char *log_file_name){
  unordered_map<string, string>::iterator it = LOG_MAP.begin(), end = LOG_MAP.end();
  unsigned int log_size = LOG_MAP.size(), count=0;
  while(it!=end){
    LOGS += "\t"+it->first +":"+it->second+"";
    ++count;
    ++it;
    if(count<log_size){
      LOGS += ", \n";
    }
  }

  string log_name(log_file_name);
  LOGS = log_name+":{\n"+LOGS+"\n}";
  ofstream foutput(log_file_name);
  foutput << LOGS;
  foutput.close();
}

// add log data to LOGS in JSON format
static void add_logs(const string header, const string log){
  LOG_MAP[header]=log;
}

// write contents to file
// If add_flag is true, this function works as the append mode.
// Otherwise, it rewrites the file.
static void write_file(const char *file, const string *contents, bool add_flag){
  ofstream f;
  
  if(!add_flag){
    f.open(file, ios::out);
  }else{
    f.open(file, ios::app);
  }
  
  f << *contents << endl;  
  f.close();
}

// fill a vector with random values by following the normal distribution
static void gen_normaldist_vector(vector<double> *vec, const double avg, const double sd){
  mt19937 mt{random_device{}()};
  normal_distribution<> dist(avg, sd);

  int size = (*vec).size();
  vector<double>::iterator it = vec->begin(), end = vec->end();
  for(int i=0; i<size; ++i){
    *(it+i) = dist(mt);
  }
}

#endif // UTIL_H
