/**
 * BNGenome.h
 * @description Bayes network genome.
 * @author David Durman, 2008
 */

#ifndef _BNGENOME_H
#define _BNGENOME_H

#include <cmath>
#include <vector>
#include <map>
#include <fstream>

#include "arrangement.h"
#include "configuration.h"

/**
 * Comment if you want to use non decomposable
 * scoring function. Actually, there is no good 
 * reason for it!
 */
#define SF_DECOMPOSABLE

/**
 * Comment if you want to mutate only one arc
 * at one step (i.e. number of mutatins will be 0 or 1). 
 * Otherwise, mutation will be performed on 
 * whole adjacency matrix.
 */
//#define MULTI_MUTATION

/**
 * Choose one scoring function.
 */
#define BIC_SCORE
//#define AIC_SCORE

//#define cout STD_COUT

/**
 * Used to shift the result of evaluation to obtain
 * positive numbers only.
 * Change if you need.
 */
#define EVALUATION_SHIFT 1000000.0

/**
 * Type of a dataset.
 */
typedef unsigned** tDataset;


class BNGenome : public GA2DBinaryStringGenome {
public:
  /**
   * Type of a probabilistic model (of it's parameters).
   *
   * configuration: 
   *	node , c(x_i, x_pa(i))
   *
   * configuration < correspondece , from_count >
   * correspondece / from_count = ML estimation of p, i.e. E( p(x_i, x_pa(i)) ) = n(x_i, x_pa(i)) / n(x_pa(i))
   *
   */
  typedef std::map<configuration, std::pair<int, int> > tProbModel;
  /**
   * Type of a model.
   * <tProbModel, int> ... probabilistic model and number of parameters of whole model.
   */
  typedef std::pair<tProbModel, int> tModel;  // probabilistic model + its size (number of parameters of BN)

public:
  GADefineIdentity("BNGenome", 201);
  static void Init(GAGenome&);
  static int Mutate(GAGenome&, float);
  /**
   * Used for diversity calculations.
   * @return 0 if two genomes are identical (no diversity), otherwise number of bits on same positions
   */
  static float Compare(const GAGenome&, const GAGenome&);
  /**
   * This is the objective function.  
   * Uses BIC/AIC scoring function to score candidates.
   */
  static float Evaluate(GAGenome&);
  static int Cross(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);

public:
  /**
   * @param int nodes number of nodes, i.e. width and hight of adjacency matrix
   * @param tDataset dataset two dimensional array of data samples
   */
  BNGenome(int nodes, tDataset dataset) : GA2DBinaryStringGenome(nodes, nodes, NULL, (void*) dataset) { 
    initializer(Init);
    comparator(Compare);
    mutator(Mutate);
    evaluator(Evaluate); 
    crossover(Cross); 
  }

  /**
   * Used to populate a population.
   */
  BNGenome(const BNGenome& orig) : GA2DBinaryStringGenome((GA2DBinaryStringGenome&) orig) { 
#ifdef SF_DECOMPOSABLE
    nodesContribs = orig.nodesContribs;
#endif
  }

  virtual ~BNGenome() { }

  /**
   * This is not used anywhere in this case.
   */
  BNGenome& operator=(const BNGenome& orig){
    if(&orig != this){
      copy(orig);
    }
    return *this;
  }

  virtual BNGenome* clone(CloneMethod) const {
    return new BNGenome(*this);
  }

  /**
   * This is not used anywhere in this case.
   */
  virtual void copy(const GAGenome& orig) {
    GA2DBinaryStringGenome::copy(orig);  // this copies all of the base genome parts
#ifdef SF_DECOMPOSABLE
    const BNGenome& o = (BNGenome&) orig;
    nodesContribs = o.nodesContribs;
#endif
  }

  virtual int equal(const GAGenome& g){
    int count = 0;
    const BNGenome &sis=(BNGenome &)*this;
    const BNGenome &bro=(BNGenome &)g;
    for (int icol = 0; icol < sis.width(); icol++)
      for (int irow = 0; irow < sis.width(); irow++)
	if (sis.gene(icol, irow) != bro.gene(icol, irow))
	  count += 1;

#ifdef SF_DECOMPOSABLE
    std::map<int, float>::const_iterator itbro;
    for (std::map<int, float>::const_iterator itsis = sis.nodesContribs.begin();
	 itsis != sis.nodesContribs.end(); itsis++){
      itbro = bro.nodesContribs.find(itsis->first);
      if (itbro == bro.nodesContribs.end() ||
	  itbro->second != itsis->second)
	count += 1;
    }
#endif

    return (count == 0);
  }


public:
  /**
   * Computes log-likelihood function evaluated at its maximum ( L(M) ) of the model.
   * @param tModel m model representation
   * @return double
   */
  double logLikelihood(tModel*);
  double logLikelihood();

  /**
   * BIC score of the model.
   */
  double BIC();
  /**
   * AIC score of the model.
   */
  double AIC();

  /**
   * Get the probability model obtained from the dataset using my BN structure.
   */
  tModel* getModel();  

public:
  //TODO: make these two attributes static!
  static int observs;	// number of observations in the dataset
  static std::vector<int> values;	// possible values for each node

#ifdef SF_DECOMPOSABLE
  // nodes and its contribution to log-likelihood
  // if a node isn't there, we have to recompute its contribution and put it in there
  std::map<int, float> nodesContribs;
#endif

  static void numberOfObservations(int o){
    observs = o;
  }
  static void possibleValuesForNodes(std::vector<int>& v){
    values = v;
  }
  /**
   * Use if all nodes have the same possible values.
   * @param int v possible values (e.g. 2 means possible values 0, 1 for each node)
   * @param int nodes number of nodes
   */
  static void possibleValuesForNodes(int v, int nodes){
    for (int i = 0; i < nodes; i++)
      values.push_back(v);
  }


  /**
   * Cycle detection part.
   * DFS with a few modifications.
   */



public:
  /**
   * @return true if my graph is a DAG
   */
  bool isDag(){ 
    cycleSearchInit();
    delete[] marked;
    delete[] inProg;
    return dagflag; 
  }  

private:  
  bool* marked;
  bool* inProg;
  bool dagflag;

private:
  // standard DFS with a few modifications
  void cycleSearchInit(){
    dagflag = true;
    marked = new bool[width()];
    inProg = new bool[width()];

    // init
    for (int i = 0; i < width(); i++){
      marked[i] = false;
      inProg[i] = false;
    }

    for (unsigned v = 0; v < (unsigned) width(); v++)
      if (!marked[v])
	cycleSearch(v);
  }

  void cycleSearch(unsigned v){
    marked[v] = true;
    inProg[v] = true;

    // for all adjacents
    std::vector<unsigned> adjs = adjacents(v);
    std::vector<unsigned>::const_iterator w = adjs.begin();
    for (; w != adjs.end(); w++){
      if (inProg[*w]){
	dagflag = false;
	return;
      } else if (!marked[*w]){
	cycleSearch(*w);
      }
    }
    inProg[v] = false;
  }
  
  // get all adjacents of node v
  std::vector<unsigned> adjacents(unsigned v){
    std::vector<unsigned> adjs;
    for (int i = 0; i < width(); i++){
      if ((gene(v, i) == 1))
	adjs.push_back(i);
    }
    return adjs;
  }

public:

  /**
   * Generates the dot file. Input for Graphviz.
   * @param string name name of the graph
   */
  void genDotFile(std::string name, bool verbose){
    if (verbose)
      std::cout << "Generating dot file " << (name + ".dot .");

    std::ofstream dotFile((name + ".dot").c_str());
    dotFile << "digraph " << name << " { " << std::endl;
    dotFile << "\tsize=\"8.5\"" << std::endl;
    dotFile << "\tnode [shape = circle]" << std::endl;

    for (int i = 0; i < width(); i++){
      for (int j = 0; j < width(); j++){
	if (verbose) std::cout << ".";
	if (gene(i, j) == 1)
	  dotFile << "\t\"node_" << i << "\" -> \"node_" << j << "\"" << std::endl;
      }
    }
    if (verbose)
      std::cout << std::endl;
    
    dotFile << "}" << std::endl;
    dotFile.close();
  }

  /**
   * Print the probabilistic model. (without already calculated contributions)
   */
  void printModel(const tModel* model, std::ostream& os){
    tProbModel::const_iterator it = model->first.begin();
    for (; it != model->first.end(); it++){
      configuration c = it->first;
      c.print(os);
      os << ":" << it->second.first << "/" << it->second.second << std::endl;
    }
    //    os << "number of parameters: " << model->second << std::endl;
  }

#ifdef SF_DECOMPOSABLE
  /**
   * Saved contribtuions to the probabilistic model
   * which are not presented in the model.
   */
  void printContributions(){
    std::cout << "contributions: ";
    std::map<int, float>::const_iterator it = nodesContribs.begin();
    for (; it != nodesContribs.end(); it++)
      std::cout << "node(" << it->first << ") = " << it->second << "  ";
    std::cout << std::endl;
  }
#endif
};

#endif
