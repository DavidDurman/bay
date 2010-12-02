/**
 * baybuild.cc
 * @description <p>Learn Bayes network structure 
 *	from a dataset using Genetic algorithms.</p>
 * @author David Durman, 2008
 */

#include <ga/GAGenome.h>
#include <ga/GA2DBinStrGenome.h>
#include <ga/std_stream.h>

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "BNDataset.h"
#include "BNGenome.h"

using namespace std;

#define cout STD_COUT

// default options
#define POPSIZE	100
#define NGEN	20
#define PMUT	0.9
#define PCROSS	0.1
#define BN_GRAPH_NAME "bn"
#define FDATASET "dataset.csv"

/**
 * Which algorithm to use
 */
//#define STEADYSTATE_GA
//#define INCREMENTAL_GA

#ifdef STEADYSTATE_GA
#include<ga/GASStateGA.h>
#elif defined(INCREMENTAL_GA)
#include<ga/GAIncGA.h>
#else
#include<ga/GASimpleGA.h>
#endif

/**
 * Options for Steady state GA.
 */
//#define PREPLACEMENT 0.2


/**
 * Options for Incremental GA.
 */
//#define REPLACEMENT GAIncrementalGA::RANDOM
//#define REPLACEMENT GAIncrementalGA::BEST
//#define REPLACEMENT GAIncrementalGA::WORST
//#define REPLACEMENT GAIncrementalGA::CUSTOM
//#define REPLACEMENT GAIncrementalGA::CROWDING
//#define REPLACEMENT GAIncrementalGA::PARENT 


/**
 * Options for GA.
 */
//#define TOURNAMENT_SELECTOR
//#define ELITIST


void printUsage(){
  cout << "Usage: baybuild [-s|-d|-p|-g|-m|-c|-h|-b|-v]" << endl;
  cout << "\t-h      \tprint program options" << endl;
  cout << "\t-s seed \tset random seed" << endl;
  cout << "\t-d file \tdataset csv file (comma separated)" << endl;
  cout << "\t-p size \tGA population size" << endl;
  cout << "\t-g count\tGA number of generations" << endl;
  cout << "\t-m prob \tGA probability of mutation" << endl;
  cout << "\t-c prob \tGA probability of crossover" << endl;
  cout << "\t-b name \tBN name; generates dot file 'name.dot'" << endl;
  cout << "\t-v      \tverbose" << endl;
}

/*--------------------------------------------------.
  |  Main.                                            |
  `--------------------------------------------------*/

int main(int argc, char **argv){

  if (argc == 1){
    printUsage();
    return 0;
  }

  /**
   * Default options.
   */

  string datasetFileName = FDATASET;
  int popsize = POPSIZE;
  int ngen = NGEN;
  float pmut = PMUT;
  float pcross = PCROSS;
  string bn_graph_name = BN_GRAPH_NAME;
  bool verbose = false;


  /**
   * Process command line options.
   */

  for(int ii = 1; ii < argc; ii++) {
    if (strcmp(argv[ii], "-s") == 0) {		// seed
      GARandomSeed((unsigned int)atoi(argv[++ii]));
    } else if (strcmp(argv[ii], "-d") == 0) {	// dataset
      datasetFileName = string(argv[++ii]);

    } else if (strcmp(argv[ii], "-p") == 0){	// population size
      popsize = atoi(argv[++ii]);
    } else if (strcmp(argv[ii], "-g") == 0){	// number of generations
      ngen = atoi(argv[++ii]);
    } else if (strcmp(argv[ii], "-m") == 0){	// probability of mutation
      pmut = atof(argv[++ii]);
    } else if (strcmp(argv[ii], "-c") == 0){	// probability of crossover
      pcross = atof(argv[++ii]);

    } else if (strcmp(argv[ii], "-b") == 0){	// BN graph name
      bn_graph_name = string(argv[++ii]);

    } else if (strcmp(argv[ii], "-v") == 0){	// be verbose
      verbose = true;

    } else if (strcmp(argv[ii], "-h") == 0){	// help
      printUsage();
      return 0;
    }
  }
  
  /**
   * Load a dataset and print its parameters if verbose flag is presented.
   */

  int nnodes, nobservs;
  datasetParams(datasetFileName, &nobservs, &nnodes);

  if (verbose){
    cout << "Number of observations: " << nobservs << endl << "Number of nodes: " << nnodes << endl;
    cout << "Population size: " << popsize << endl;
    cout << "Number of generations: " << ngen << endl;
    cout << "Probability of mutation: " << pmut << endl;
    cout << "Probability of crossover: " << pcross << endl;
    cout << "Graph name: " << bn_graph_name << endl;
    cout << "Dataset: " << datasetFileName << endl;
    cout << "--------------------------------" << endl;
    cout << endl;
  }

  tDataset dataset = loadDataset(datasetFileName, nobservs, nnodes);
  // possible values for each node
  vector<int> values = datasetNodeValues(dataset, nobservs, nnodes);

  /**
   * Prepare and run genetic algorithm.
   */

  BNGenome genome(nnodes, dataset);
  genome.numberOfObservations(nobservs);
  genome.possibleValuesForNodes(values);

#ifdef STEADYSTATE_GA
  GASteadyStateGA ga(genome);
  ga.pReplacement( PREPLACEMENT );
#elif defined(INCREMENTAL_GA)
  GAIncrementalGA ga(genome);
  ga.replacement( REPLACEMENT );
#else
  GASimpleGA ga(genome);  
#endif

  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);

#ifdef TOURNAMENT_SELECTOR
  GATournamentSelector* sel = new GATournamentSelector;
  ga.selector(*sel);
#endif

#ifdef ELITIST
  ga.elitist(gaTrue);
#endif

  ga.initialize();
  ga.evolve(); 

#ifdef TOURNAMENT_SELECTOR
  delete sel;	
#endif


  /**
   * Get solution and print results if verbose flag is presented.
   */

  BNGenome& theBestOne = (BNGenome &) ga.statistics().bestIndividual();
  // probability model
  theBestOne.nodesContribs.clear();	// to obtain the entire probability model
  BNGenome::tModel* pm = theBestOne.getModel();

  if (verbose){
    // print the results
    cout << "Number of parameters of the best model: " << pm->second << endl;
    cout << "Log-likelihood of the best model: ";
    double ll = theBestOne.logLikelihood(pm);
    cout <<  ll << endl;

    cout << "BIC score of the best model: ";
    cout << ll - (log(theBestOne.observs) / 2) * pm->second << endl << endl;

    cout << "Structure of the best model: " << endl;
    cout << theBestOne;
    cout << endl << endl;
    
    cout << "Probability model:" << endl;
    theBestOne.printModel(pm, cout);
    cout << endl;
  }

  // generate "dot" file for graphviz
  theBestOne.genDotFile(bn_graph_name, verbose);

  // save topological structure of the best BN
  cout << "Saving topological structure of the best BN to file " << bn_graph_name + ".net" << endl;
  std::ofstream topBNfile((bn_graph_name + ".net").c_str());
  topBNfile << theBestOne;
  topBNfile.close();

  // save probability model of the best BN
  cout << "Saving probability model of the best BN to file " << bn_graph_name + ".prob" << endl;
  std::ofstream probBNfile((bn_graph_name + ".prob").c_str());
  theBestOne.printModel(pm, probBNfile);
  probBNfile.close();

  /**
   * Clean the rubish.
   */

  BNGenome::tProbModel::iterator iter = pm->first.begin();
  for (; iter != pm->first.end(); iter++)
    delete [] iter->first.buff;
  delete pm;
  freeDataset(dataset, nobservs);

  return 0;
}
  

