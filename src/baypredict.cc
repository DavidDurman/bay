/**
 * baypredict.cc
 * @description <p>Builds BN from topological
 *	structure information and probability model
 *	and predict diagnose from evidences using
 *	BN inference</p>
 * @author David Durman, 2008
 */

#include "dlib/bayes_utils.h"
#include "dlib/graph_utils.h"
#include "dlib/graph.h"
#include "dlib/directed_graph.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <set>
#include <iterator>
#include <algorithm>
#include <cassert>


using namespace dlib;
using namespace bayes_node_utils;

// default values
#define FBNTOPOLOGY "bn.net"		// file with BN adjacency matrix
#define FBNPROBABILITIES "bn.prob"	// file with BN probability model

// First we need to create an undirected graph which contains set objects at each node and
// edge.  This long declaration does the trick.
typedef graph<dlib::set<unsigned long>::compare_1b_c, dlib::set<unsigned long>::compare_1b_c>::kernel_1a_c join_tree_type;

/**
 * Probability model.
 * <node, pair<configuration, probability> >
 * configuration ... [myvalue, parents_values...]
 */
typedef std::vector<int> tConfiguration;
typedef std::pair<tConfiguration, float>  tAssessement;
typedef std::multimap<int, tAssessement> tProbModel;

/**
 * BN graph utility to find connected components and
 * their nodes, etc...
 */
class BNGraphUtility {
private:
  int components;
  int* cc;
  int** am;	// adjacency matrix
  int nnodes;	// number of nodes

public:
  BNGraphUtility(int** m, int n){
    am = m;
    nnodes = n;
    components = 1;
    cc = new int[nnodes];
    for (int i = 0; i < nnodes; i++)
      cc[i] = -1;	// unmarked
    for (int i = 0; i < nnodes; i++)
      if (cc[i] == -1){
	dfs(i);
	components++;
      }
  }

  void dfs(int node){
    cc[node] = components;
    std::vector<unsigned> adjs = adjacents(node);
    std::vector<unsigned>::const_iterator it = adjs.begin();
    for (; it != adjs.end(); it++){
      if (cc[*it] == -1){
	dfs(*it);
      }
    }
  }

  /**
   * Is there an undirected path between node v and w?
   */
  bool connected(unsigned v, unsigned w){
    return cc[v] == cc[w];
  }

  /**
   * Get number of components.
   */  
  int nComponents(){ return components; }

  /**
   * Get number of nodes in component c.
   */
  int nNodesInComponent(int c){
    int ret = 0;
    for (int i = 0; i < nnodes; i++)
      if (cc[i] == c) ret++;
    return ret;
  }

  /**
   * Get all nodes from component c.
   */
  std::vector<int> nodesInComponent(int c){
    std::vector<int> nodes;
    for (int i = 0; i < nnodes; i++)
      if (cc[i] == c) 
	nodes.push_back(i);
    return nodes;
  }

  /**
   * To which component does node v belong?
   */
  int componentForNode(int v){
    return cc[v];
  }

  /**
   * Get all adjacents of node v.
   */
  std::vector<unsigned> adjacents(unsigned v){
    std::vector<unsigned> adjs;

    // find adjacents
    for (int i = 0; i < nnodes; i++){
      if ((am[i][v] == 1) || (am[v][i] == 1))
	adjs.push_back(i);
    }
    return adjs;
  }

  ~BNGraphUtility(){ delete [] cc; }
};


/**
 * Print probability model. (Just for debugging purposes)
 */
void printModel(const tProbModel& m){
  tProbModel::const_iterator it = m.begin();
  for (; it != m.end(); it++){
    int n = it->first;
    const std::vector<int>& c = it->second.first;
    float p = it->second.second;

    std::cout << n << ":";
    std::vector<int>::const_iterator cit = c.begin();
    for (; cit != c.end(); cit++)
      std::cout << *cit << " ";
    std::cout << ":" << p << std::endl;
  }
}

void cleanAdjacencyMatrix(int** am, int nnodes){
  // deallocate adjacency matrix
  for (int i = 0; i < nnodes; i++)
    delete [] am[i];
  delete [] am;
}

void printUsage(){
  std::cout << "Usage: baypredict [-n|-p|-h]" << std::endl;
  std::cout << "\t-h      \tprint program options" << std::endl;
  std::cout << "\t-n file \tBN adjacency matrix file" << std::endl;
  std::cout << "\t-p file \tBN probability model file" << std::endl;      
  std::cout << "\t-e evidences \tcomma separated evidences assignments (e.g. 0=1,5=1,3=0 means\n\t\t\t evidence node 0 is equal to 1 (true) is observed, likewise for the others)" << std::endl;      
  std::cout << "\t-c cause \tcause assignment (e.g. 0=1 means what is the probability of node 0 being true?" << std::endl;            
}

/**
 * Main.
 */

int main(int argc, char** argv){

  if (argc == 1){
    printUsage();
    return 0;
  }

  std::string topBNfilename(FBNTOPOLOGY);
  std::string probBNfilename(FBNPROBABILITIES);

  // <node, value>
  std::map<int, int> evidences;	// evidences
  int causeNode = 3;		// cause node
  int causeValue = 1;		// cause value


  /**
   * Process command line options.
   */
  
  for (int ii = 1; ii < argc; ii++){
    if (strcmp(argv[ii], "-n") == 0){
      topBNfilename = std::string(argv[++ii]);

    } else if (strcmp(argv[ii], "-p") == 0){
      probBNfilename = std::string(argv[++ii]);


      /**
       * Parsing evidences.
       */
    } else if (strcmp(argv[ii], "-e") == 0){
      std::string sEvidences = std::string(argv[++ii]);
      std::istringstream sEvidencesStream(sEvidences);
      std::string sEvidence;
      while (getline(sEvidencesStream, sEvidence, ',')){
	std::istringstream sEvidenceStream(sEvidence);
	std::string sEvidenceNode;
	std::string sEvidenceValue;
	if (getline(sEvidenceStream, sEvidenceNode, '=')){
	  if (getline(sEvidenceStream, sEvidenceValue, '=')){
	    std::istringstream sEvidenceNodeStream(sEvidenceNode);
	    std::istringstream sEvidenceValueStream(sEvidenceValue);
	    int evidenceNode;
	    int evidenceValue;
	    sEvidenceNodeStream >> evidenceNode;
	    sEvidenceValueStream >> evidenceValue;
	    evidences[evidenceNode] = evidenceValue;
	  } else {
	    std::cerr << "Error in parsing evidence value." << std::endl;	    
	    return 1;
	  }
	} else {
	  std::cerr << "Error in parsing evidence node." << std::endl;
	  return 1;
	}
      }
      
      /**
       * Parsing cause.
       */
    } else if (strcmp(argv[ii], "-c") == 0){
      std::string sCause = std::string(argv[++ii]);
      std::istringstream sCauseStream(sCause);
      std::string sCauseNode;
      std::string sCauseValue;
      if (getline(sCauseStream, sCauseNode, '=')){
	if (getline(sCauseStream, sCauseValue, '=')){
	  std::istringstream sCauseNodeStream(sCauseNode);
	  std::istringstream sCauseValueStream(sCauseValue);
	  sCauseNodeStream >> causeNode;
	  sCauseValueStream >> causeValue;
	} else {
	  std::cerr << "Error in parsing cause node." << std::endl;
	  return 1;
	}
      } else {
	std::cerr << "Error in parsing cause node." << std::endl;
	return 1;
      }
      
    } else if (strcmp(argv[ii], "-h") == 0){
      printUsage();
      return 0;
    }
  }

  /*
    std::cout << "Evidences:" << std::endl;
    std::map<int, int>::const_iterator eit = evidences.begin();
    for (; eit != evidences.end(); eit++)
    std::cout << "\tnode " << eit->first << " = " << eit->second << " is observed" << std::endl;
    std::cout << "Given these evidences, what is the probability of node " << causeNode << " = " << causeValue << " ?" << std::endl;
  */

  int nnodes = 0;	// number of nodes

  /**
   * Load probability model.
   */

  tProbModel pModel;	// probabilistic model
  std::map<int, std::set<int> > nodeNumValues;

  std::ifstream probBNfile(probBNfilename.c_str());
  std::string line;
  while (getline(probBNfile, line)){
    std::istringstream linestream(line);
    std::string sNode;
    std::string sConfiguration;
    std::string sProbability;

    getline(linestream, sNode, ':');
    getline(linestream, sConfiguration, ':');
    getline(linestream, sProbability, ':');

    // node
    int node;
    std::istringstream nodeStream(sNode);
    nodeStream >> node;

    // configuration
    std::vector<int> configuration;
    std::string sConfigurationItem;
    std::istringstream configurationStream(sConfiguration);

    // first my value (always present)
    getline(configurationStream, sConfigurationItem, ' ');
    std::istringstream configurationItemStream(sConfigurationItem);
    int configurationItem;
    configurationItemStream >> configurationItem;
    nodeNumValues[node].insert(configurationItem);
    configuration.push_back(configurationItem);
    
    // parents values
    while (getline(configurationStream, sConfigurationItem, ' ')){
      std::istringstream configurationItemNextStream(sConfigurationItem);
      configurationItemNextStream >> configurationItem;
      configuration.push_back(configurationItem);
    }

    // probability
    unsigned pMe, pAll;
    std::istringstream probabilityStream(sProbability);
    std::string spMe;
    std::string spAll;
    getline(probabilityStream, spMe, '/');
    getline(probabilityStream, spAll, '/');

    std::istringstream spMeStream(spMe);
    std::istringstream spAllStream(spAll);
    spMeStream >> pMe;
    spAllStream >> pAll;    

    float probability = (float) pMe / pAll;

    // increase number of nodes if this one is not already in the probability model
    if (pModel.find(node) == pModel.end())
      nnodes++;

    pModel.insert(std::pair<
		  int, 
		  std::pair<std::vector<int>, float> 
		  >
		  (
		   node, 
		   std::pair<std::vector<int>, float>(configuration, probability)
		   )
		  );

  }
  probBNfile.close();

  //  printModel(pModel);

  /**
   * Load topological structure of BN.
   */

  // <node, [its parents...]>
  std::map<int, std::vector<int> > parents;

  // allocate adjacency matrix
  int** am = new int*[nnodes];
  for (int i = 0; i < nnodes; i++){
    am[i] = new int[nnodes];
  }

  unsigned linenum = 0;
  std::ifstream topBNfile(topBNfilename.c_str());
  
  while (getline(topBNfile, line)){
    std::istringstream linestream(line);
    std::string item;

    std::vector<char> l;
    copy(std::istream_iterator<char>(linestream),
	 std::istream_iterator<char>(), back_inserter(l));
    
    // for all nodes
    for (int i = 0; i < l.size(); i++){
      // add an edge if there is any
      if (l[i] == '1'){
	am[i][linenum] = 1;
	parents[linenum].push_back(i);
      } else {
	am[i][linenum] = 0;
      }
    }
    linenum++;
  }
  topBNfile.close();


  /**
   * Prediction.
   */

  BNGraphUtility bnUtil(am, nnodes);

  // consider nodes from the same component only
  int causeComponent = bnUtil.componentForNode(causeNode);
  std::vector<int> causeComponentNodes = bnUtil.nodesInComponent(causeComponent);

  /**
   * If there is only one node, just print its probability 
   * and terminate.
   */

  if (causeComponentNodes.size() == 1){
    tProbModel::const_iterator pmit = pModel.find(causeNode);
    assert(pmit != pModel.end());
    
    do {
      tAssessement assessment = pmit->second;
      tConfiguration configuration = assessment.first;
      float probability = assessment.second;
      int nodeVal = configuration[0];	// always presented
      
      if (nodeVal == causeValue){
	std::cout << "p(" << causeNode << " = " << causeValue << ") = " << probability << std::endl;	
	cleanAdjacencyMatrix(am, nnodes);
	return 0;
      }

      pmit++;
    } while (pmit != pModel.upper_bound(causeNode));
  }

  // bayes network
  directed_graph<bayes_node>::kernel_1a_c bn;
  bn.set_number_of_nodes(nnodes);

  /**
   * For all nodes in this component.
   */

  std::vector<int>::const_iterator node = causeComponentNodes.begin();
  for (; node != causeComponentNodes.end(); node++){

    /**
     * Set edges.
     */

    std::vector<int>::const_iterator pit = parents[*node].begin();
    // for all parents of this node
    for (; pit != parents[*node].end(); pit++){
      bn.add_edge(*pit, *node);
    }

    /**
     * Set number of values for each node.
     */
      
    set_node_num_values(bn, *node, nodeNumValues[*node].size());


    /**
     * Assign probabilities to all nodes.
     */

    // for all nodes
    assignment parent_state;

    tProbModel::const_iterator pmit = pModel.find(*node);
    if (pmit != pModel.end())
      do {
	tAssessement assessment = pmit->second;
	tConfiguration configuration = assessment.first;
	float probability = assessment.second;
	int nodeVal = configuration[0];	// always presented

	parent_state.clear();

	tConfiguration::const_iterator cit = configuration.begin();
	int icit = 0;	// position in configuration (without node)
	// set parent state according to the configuration
	if (parents.find(*node) != parents.end()){
	  for (++cit; cit != configuration.end(); cit++, icit++){
	    parent_state.add(parents[*node].at(icit), *cit);
	  }
	}
	set_node_probability(bn, *node, nodeVal, parent_state, probability);

	pmit++;
      } while (pmit != pModel.upper_bound(*node));
  }//endfor all nodes in this component

  join_tree_type join_tree;

  // Now we need populate the join_tree with data from our bayesian network.  The next to 
  // function calls do this.  Explaining exactly what they do is outside the scope of this
  // example.  Just think of them as filling join_tree with information that is useful 
  // later on for dealing with our bayesian network. 
  create_moral_graph(bn, join_tree);

  create_join_tree(join_tree, join_tree);
  
  // Now we have a proper join_tree we can use it to obtain a solution to our
  // bayesian network.  Doing this is as simple as declaring an instance of
  // the bayesian_network_join_tree object as follows:
  bayesian_network_join_tree solution(bn, join_tree);


  /**
   * Print the result.
   */

  if (evidences.size() == 0){
    // no evidences used
    std::cout << "p(" << causeNode << " = " << causeValue << ") = " << solution.probability(causeNode)(causeValue) << std::endl;

  } else {
    
    std::cout << "p(" << causeNode << " = " << causeValue;

    bool eFirst = true;	// just for fancy print
    std::map<int, int>::const_iterator evid = evidences.begin();
    for (; evid != evidences.end(); evid++){
      if (!bnUtil.connected(evid->first, causeNode))
	continue;

      // tell BN about evidences
      set_node_value(bn, evid->first, evid->second);
      set_node_as_evidence(bn, evid->first);

      if (!eFirst)
	std::cout << ", ";
      else {
	std::cout << " | ";	
	eFirst = false;
      }

      std::cout << evid->first << " = " << evid->second;
    }

    bayesian_network_join_tree solution_with_evidence(bn, join_tree);
    std::cout << ") = " << solution_with_evidence.probability(causeNode)(causeValue) << std::endl;

  }

  /**
   * Clean the rubish.
   */

  cleanAdjacencyMatrix(am, nnodes);

  return 0;
}
