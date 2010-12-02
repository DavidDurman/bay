/**
 * BNDataset.h
 * @description Loading dataset and its parameters.
 * @author David Durman, 2008
 */

#ifndef _BNDATASET_
#define _BNDATASET_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

/*--------------------------------------------------.
|  Data type                                        |
`--------------------------------------------------*/

/**
 * Type of a dataset.
 */
typedef unsigned** tDataset;

/*--------------------------------------------------.
|  Prototypes                                       |
`--------------------------------------------------*/

/**
 * Get number of observations and number of nodes from csv.
 */
void datasetParams(string filename, int* nobservs, int* nodes);
/**
 * Get possible values for each node.
 */
vector<int> datasetNodeValues(tDataset dataset, int nobservs, int nodes);
/** 
 * Loads dataset from csv file.
 * @param string filename path to a csv file
 * @param unsigned observations number of rows
 * @param unsigned nodes number of columns
 * @return tDataset dataset as 2D array of unsigned integers
 */
tDataset loadDataset(string filename, unsigned observations, unsigned nodes);
/**
 * Free dataset memory.
 */
void freeDataset(tDataset dataset, unsigned observations);


/*--------------------------------------------------.
|  Implementation                                   |
`--------------------------------------------------*/

/**
 * Get number of observations and number of nodes from csv.
 */
void datasetParams(string filename, int* nobservs, int* nodes){
  *nobservs = 0;
  *nodes = 0;

  ifstream file(filename.c_str());
  string line;

  // first line -> get number of nodes
  if (!getline(file, line))
    return;

  (*nobservs)++;

  istringstream linestream(line);
  string item;
  while (getline(linestream, item, ','))
    (*nodes)++;

  // remaining lines
  while (getline(file, line)){
    (*nobservs)++;
  }
}

/**
 * Get possible values for each node.
 */
vector<int> datasetNodeValues(tDataset dataset, int nobservs, int nodes){
  vector< set<int> > vals(nodes, set<int>() );
  for (int n = 0; n < nodes; n++){
    for (int o = 0; o < nobservs; o++){
      vals[n].insert(dataset[o][n]);
    }
  }

  vector<int> values;
  vector< set<int> >::const_iterator iter = vals.begin();
  for (; iter != vals.end(); iter++)
    values.push_back(iter->size());

  return values;
}

/** 
 * Loads dataset from csv file.
 * @param string filename path to a csv file
 * @param unsigned observations number of rows
 * @param unsigned nodes number of columns
 * @return tDataset dataset as 2D array of unsigned integers
 */
tDataset loadDataset(string filename, unsigned observations, unsigned nodes){
  // allocate dataset memory
  tDataset dataset = new unsigned*[observations];
  for (unsigned i = 0; i < observations; i++)
    dataset[i] = new unsigned[nodes];

  ifstream file(filename.c_str());
  string line;
  unsigned linenum = 0;
  
  while (getline(file, line)){
    istringstream linestream(line);
    string item;
    int itemnum = 0;
    
    while (getline(linestream, item, ',')){
      // convert string value into integer value
      istringstream auxstream(item);
      int value;
      auxstream >> value;
      // assign value to the dataset
      dataset[linenum][itemnum] = value;
      itemnum++;
    }
    linenum++;
  }
  return dataset;
}

/**
 * Free dataset memory.
 */
void freeDataset(tDataset dataset, unsigned observations){
  for (unsigned i = 0; i < observations; i++)
    delete [] dataset[i];
  delete [] dataset;
}



#endif
