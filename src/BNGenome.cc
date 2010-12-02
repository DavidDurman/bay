/**
 * BNGenome.cc
 * @description Bayes network genome.
 * @author David Durman, 2008
 */

#include <ga/GAGenome.h>
#include <ga/GA2DBinStrGenome.h>
#include <ga/std_stream.h>

#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "BNGenome.h"

using namespace std;

/**
 * Used to init first population.
 */
void BNGenome::Init(GAGenome& c){
  //  cout << "init" << endl;
  BNGenome& child = (BNGenome &) c;

  for (int icol = 0; icol < child.width(); icol++){
    for (int irow = 0; irow < child.width(); irow++)    
      child.gene(icol, irow, 0);
  }

  int icol = GARandomInt(0, child.width() - 1);
  int irow = GARandomInt(0, child.width() - 1);

  // put one arc with some probability
  if (GAFlipCoin(1.0) && (icol != irow))
    child.gene(icol, irow, 1);

  child._evaluated = gaFalse;

#ifdef SF_DECOMPOSABLE
  child.nodesContribs.clear();
#endif
}

/**
 * Mutation operator.
 * Add/remove an arc/arcs.
 */
int BNGenome::Mutate(GAGenome& c, float p){
  //  cout << "mutate" << endl;
  BNGenome& child = (BNGenome &) c;  
  bool change = false;		// any change happend?
  int nMut = 0;			// number of mutations
  child._evaluated = gaTrue;    // suppose no change

#ifdef MULTI_MUTATION
  for (int icol = 0; icol < child.width(); icol++)
    for (int irow = 0; irow < child.width(); irow++)
#else
  int icol = GARandomInt(0, child.width() - 1);
  int irow = GARandomInt(0, child.width() - 1);
#endif

  if (GAFlipCoin(p) && (icol != irow)){
    child.gene(icol, irow, ((child.gene(icol, irow) == 0) ? 1 : 0));
    change = true;

    // prevent of bidirectional arcs and cyclic graphs
    if ( ((child.gene(icol, irow) == child.gene(irow, icol)) && (child.gene(icol, irow) == 1)) ||
	 !child.isDag() ){
      child.gene(icol, irow, ((child.gene(icol, irow) == 0) ? 1 : 0));
      change = false;
    }

    if (change){
#ifdef SF_DECOMPOSABLE
      child.nodesContribs.erase(irow);
#endif
      nMut++;
      child._evaluated = gaFalse;
    } 
  }

  assert(child.isDag() == true);  
  return nMut;
}

/**
 * Crossover operator. 
 * Randomly choose one node and swap columns 
 * of that node across parents in adjacency matrix.
 * This leads to two offsprings.
 */
int BNGenome::Cross(const GAGenome& p1, const GAGenome& p2,
		    GAGenome* c1, GAGenome* c2){
  //  cout << "cross" << endl;
  BNGenome& mom = (BNGenome &) p1;
  BNGenome& dad = (BNGenome &) p2;

  int nCross = 0;
  int node = GARandomInt(0, mom.width() - 1);  
  
  if (c1){
    BNGenome& sis = (BNGenome &) *c1;
    sis.copy(mom);
    sis._evaluated = gaFalse;

    for (int irow = 0; irow < sis.width(); irow++){
      // change parents?
      if (sis.gene(node, irow) != dad.gene(node, irow)){
	sis.gene(node, irow, dad.gene(node, irow));
#ifdef SF_DECOMPOSABLE
	sis.nodesContribs.erase(irow);
#endif
      }
    }
    nCross++;

    // if not DAG then get back mom's graph
    if (sis.isDag() == false){
      sis.copy(mom);
      nCross--;
      sis._evaluated = gaTrue;
    }
  }

  if (c2){
    BNGenome& bro = (BNGenome &) *c2;    
    bro.copy(dad);
    bro._evaluated = gaFalse;

    for (int irow = 0; irow < bro.width(); irow++){
      // change parents?
      if (bro.gene(node, irow) != mom.gene(node, irow)){
	bro.gene(node, irow, mom.gene(node, irow));
#ifdef SF_DECOMPOSABLE
	bro.nodesContribs.erase(irow);
#endif
      }
    }
    nCross++;

    // if not DAG then get back mom's graph
    if (bro.isDag() == false){
      bro.copy(dad);
      nCross--;
      bro._evaluated = gaTrue;
    }
  }
  return nCross;
}



/**
 * Used for computing diversity.
 * Number of non-equally bits between sis and bro.
 */
float BNGenome::Compare(const GAGenome& a, const GAGenome& b){
  float count = 0.0;
  const BNGenome &sis=(BNGenome &)a;
  const BNGenome &bro=(BNGenome &)b;
  for (int icol = 0; icol < sis.width(); icol++)
    for (int irow = 0; irow < sis.width(); irow++)
      if (sis.gene(icol, irow) != bro.gene(icol, irow))
	count += 1.0;

#ifdef SF_DECOMPOSABLE
  map<int, float>::const_iterator itbro;
  for (map<int, float>::const_iterator itsis = sis.nodesContribs.begin();
       itsis != sis.nodesContribs.end(); itsis++){
    itbro = bro.nodesContribs.find(itsis->first);
    if (itbro == bro.nodesContribs.end() ||
	itbro->second != itsis->second)
      count += 1.0;
  }
#endif

  return count;
}

/**
 * Objective function. Note that it is not a fitness function.
 * Uses BIC / AIC scores.
 */
float BNGenome::Evaluate(GAGenome& g){
  //  cout << "evaluate" << endl;
  BNGenome& genome = (BNGenome &) g;  

#ifdef BIC_SCORE
  double score = genome.BIC();
#endif
#ifdef AIC_SCORE
  double score = genome.AIC();
#endif

  return score + EVALUATION_SHIFT;	
}


/**
 * Get a probabilistic model of my graph.
 * Used to compute Log-likelihood.
 */
BNGenome::tModel* BNGenome::getModel(){
  assert(isDag() == true);

  tDataset dataset = (tDataset) userData();

  map<int, vector<int> > parents;
  vector<int>::const_iterator iPar;
  tModel* model = new tModel;	// probabilistic model
  model->second = 0;		// number of parameters

  for (int icol = 0; icol < width(); icol++){
    for (int irow = 0; irow < width(); irow++){
      if (gene(icol, irow) == 1)
	parents[irow].push_back(icol); // icol is a parent of irow
    }
  }

  // for all nodes
  for (int node = 0; node < width(); node++){
    int nParents = parents[node].size();
    // how many parameters for this node in BN
    // product of all possible values for that node and all its parents
    int nParameters = 1;	


    // Number of possible values for me and all my parents in order.
    unsigned* nPossibleValuesForNode = new unsigned[1 + nParents];
    nPossibleValuesForNode[0] = values[node];

    // compute number of local parameters at the same time
    nParameters *= nPossibleValuesForNode[0];
    for (int i = 1; i < nParents + 1; i++){
      nPossibleValuesForNode[i] = values[parents[node].at(i - 1)];
      nParameters *= nPossibleValuesForNode[i];
    }
    model->second += nParameters;	// number of local parameters

#ifdef SF_DECOMPOSABLE
    // skip nodes for which parents haven't changed
    // their contribution is already calculated
    // -> don't put them to the model
    if (nodesContribs.find(node) != nodesContribs.end()){
      delete [] nPossibleValuesForNode;
      continue;
    }
#endif


    // make all possible combinations parent nodes values
    // suppose small number of parents  (otherwise combinatorial explosion!!!)
    Arrangement* arrangements = new Arrangement(nPossibleValuesForNode, nParents + 1);

    // is the arrangement wrt my parents contained in a row?
    bool correspond;	

    /**
     * For all possible combinations of parents values.
     */

    do { unsigned* a = arrangements->get();

    /**
     * For all observations.
     */

    for (int i = 0; i < observs; i++){
      unsigned* row = dataset[i];

      // Corresponds the node's value with its value in this row?
      // - if yes, check also all my parents
      correspond = true;
      int apar = 1;
      for (iPar = parents[node].begin(); iPar != parents[node].end(); ++iPar, ++apar){
	if (row[*iPar] != a[apar]){
	  correspond = false;
	  break;
	}
      }

      // if correspond, number of observations wrt my parents is incremented
      if (correspond){
	// save the configuration
	configuration c(a, 1 + nParents, node);

	// if configuration's already there, deallocate configuration's memory
	bool cDealloc = false;
	if (model->first.find(c) != model->first.end())
	  cDealloc = true;
	  
	model->first[c].second++;	// number of observations where parents have same value as the arrangement is is incremented	    

	// Even my value corresponds. Hence, number of observations wrt my parents AND me is incremented.
	if (row[node] == a[0])
	  model->first[c].first++;

	// dealloc configuration's memory
	if (cDealloc)
	  c.clear();

      }//endif correspond
    }//endfor all observations

    } while (arrangements->next());

    delete arrangements;
    delete [] nPossibleValuesForNode;

  }//end for all nodes

  return model;
}

/**
 * Computes log-likelihood function evaluated at its maximum ( L(M) ) of the model.
 * @param tModel m model representation
 * @return double
 */
double BNGenome::logLikelihood(tModel* m){
  double ll = 0.0;
  int nodes = 0;

#ifdef SF_DECOMPOSABLE
  map<int, float>::const_iterator cit;
  map<int, float> contributions;
#endif

  // compute contributions from the probabilistic model
  tProbModel::iterator fit;
  for (fit = m->first.begin(); fit != m->first.end(); fit++){

#ifdef SF_DECOMPOSABLE
    // a node can not be in my contributions and in model at the same time
    assert(nodesContribs.find(fit->first.node) == nodesContribs.end());
#endif

    double nme = (double) fit->second.first;		// n(x_i, x_pa(i))
    double nall = (double) fit->second.second;	// n(x_pa(i))
    double p = (nall == 0) ? 0.0 : (nme / nall);
    float contribution = 0.0;
    if (nme > 0.0 && p > 0.0)			// avoiding log(0) = -inf  -> use convention log(0) = 0
      contribution = nme * log(p);		// using natural logarithm

#ifdef SF_DECOMPOSABLE
    // new calculated contributions (for later usage)
    if (contributions.find(fit->first.node) != contributions.end())
      contributions[fit->first.node] += contribution;	
    else {
      contributions[fit->first.node] = contribution;
      nodes++;
    }
#endif

    ll += contribution;
  }  


#ifdef SF_DECOMPOSABLE
  // here are the remaining contributions which weren't in the model
  for (cit = nodesContribs.begin(); cit != nodesContribs.end(); cit++)
    ll += cit->second;
  assert((nodes + nodesContribs.size()) == values.size());

  // save new calculated contributions
  for (cit = contributions.begin(); cit != contributions.end(); cit++)
    nodesContribs.insert(pair<int, float>(cit->first, cit->second));
#endif

  //  cout << *this << "---" << endl;
  //  cout << "prob model:" << endl;
  //  printModel(m);
  //  printContributions();
  //  cout << "log-likelihood: " << ll << endl << endl;
  return ll;
}


double BNGenome::logLikelihood(){
  tModel* m = getModel();

  double ll = logLikelihood(m);

  // deallocate model, we no longer need it
  tProbModel::iterator it = m->first.begin();
  for (; it != m->first.end(); it++)
    delete [] it->first.buff;
  delete m;
    
  return ll;
}

double BNGenome::BIC(){
  tModel* m = getModel();

  double ll = logLikelihood(m);

  // deallocate model, we no longer need it
  tProbModel::iterator it = m->first.begin();
  for (; it != m->first.end(); it++)
    delete [] it->first.buff;
  delete m;
    
  return ll - (log(observs)/2) * m->second;
}


double BNGenome::AIC(){
  tModel* m = getModel();

  double ll = logLikelihood(m);

  // deallocate model, we no longer need it
  tProbModel::iterator it = m->first.begin();
  for (; it != m->first.end(); it++)
    delete [] it->first.buff;
  delete m;
    
  return ll - m->second;
}

// initialization of static attributes
int BNGenome::observs = 0;
vector<int> BNGenome::values = vector<int>();
