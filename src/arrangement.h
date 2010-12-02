/**
 * arrangement.h
 * @description <p>Computes all possible combinations
 *	 of node's parents values.</p>
 * @author David Durman, 2007
 */

#ifndef _ARRANGEMENT_H
#define _ARRANGEMENT_H

#include <iostream>

//////////////////
// Arrangement
class Arrangement {
public:
  unsigned *buff;
  unsigned* total_size;
  unsigned arrange_size;

  void allocate(unsigned* values, unsigned k);
  void allocate(unsigned values, unsigned k);
public:

  unsigned* limit(){ return total_size; }
  unsigned size(){ return arrange_size; }

  unsigned* get(){ return buff; }
   
  Arrangement &first(); /// < set as first permutation (identity)
  Arrangement &last(); /// < set as last permutation
  bool isfirst();
  bool islast();
  bool next(); /// < goto next permutation
  bool preview(); /// < goto preview permutation

  /**
   * @param array values possible values for all nodes
   */   
  Arrangement(unsigned* values, unsigned k = 0){ 
    //    buff = NULL;
    //    total_size = NULL;
    //    arrange_size = 0;
    allocate(values, k); first();
  }
  /**
   * @param unsigned val all nodes have the same number of values
   * TODO: allocation can not be here, i don't know when to deallocate. (see above)
   */
  Arrangement(unsigned val, unsigned k = 0){ 
    unsigned* values = new unsigned[k];
    for (unsigned i = 0; i < k; i++)
      values[i] = val;
    allocate(values, k); first();
  }
  ~Arrangement() {
    if(arrange_size){
      delete [] buff;
      //      delete [] total_size;
    }
  }
   
  Arrangement & operator=(Arrangement const &a);
  unsigned operator[](unsigned i){ 
    if(i >= arrange_size) 
      std::cerr << "operator[] ... out of range" << std::endl;
    return buff[i];
  }

  bool operator==(Arrangement const &a);
  bool operator<(Arrangement const &a);
  Arrangement &operator++(){
    if(!next()) 
      std::cerr << "operator++ ... srleady last element" << std::endl;
    return *this;
  }

  Arrangement &operator--(){
    if(!preview()) 
      std::cerr << "operator++ ... srleady first element" << std::endl;
    return *this;
  }		
};//end of class Arrangement


#endif // _ARRANGEMENT_H
