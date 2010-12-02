/**
 * arrangement.cc
 * @description <p>Computes all possible combinations
 *	 of node's parents values.</p>
 * @author David Durman, 2007
 */

#include "arrangement.h"
#include <cstring>

void Arrangement::allocate(unsigned* values, unsigned k){
  //  total_size = arrange_size = 0;
  //  if(k > n) error("allocate()", "k > n");
  if(k == 0) return;
  buff = new unsigned[k];
  if(buff) {
    total_size = values;
    arrange_size = k;
  } else {
    std::cerr << "allocate() ... Can't allocate memory" << std::endl;
  }
}

// set as First arrangement
Arrangement &Arrangement::first(){
  for(unsigned i = 0; i < arrange_size; i++)
    buff[i] = 0;
  return *this;
}

// set to last arrangement
Arrangement &Arrangement::last(){
  for(unsigned i = 0; i < arrange_size; i++)
    buff[i] = total_size[i] - 1;
  return *this;
}

bool Arrangement::isfirst(){
  for(unsigned i = 0; i < arrange_size; i++)
    if(buff[i] != 0)
      return false;
  return true;
}

bool Arrangement::islast(){
  for(unsigned i = 0; i < arrange_size; i++)
    if(buff[i] != total_size[i] - 1)
      return false;
  return true;
}

// Make next arrangement
// Return true if sucess (if not last)
bool Arrangement::next(){
  unsigned i;

  // simply increase watch
  
  i = arrange_size - 1;
  // search for incrementable marker
  while ((i > 0) && (buff[i] == total_size[i] - 1))
    i--;
  if((i == 0) && (buff[0] == total_size[i] - 1)) // not found
    return false;
  else {
    buff[i]++; //  increment
    i++;	 
    // memset(&buff[i], 0, (arrange_size-i)*sizeof(unsigned));
    for(; i < arrange_size; i++)
      buff[i] = 0;
    return true; 
  }
}

// Make preview combinations 
// Return true on sucess
bool Arrangement::preview(){
  unsigned i;

  // simply decrement watch
  
  i = arrange_size - 1;
  // search for decreasable marker
  while((i > 0) && (buff[i] == 0))
    i--;
  if((i == 0) && (buff[0] == 0))
    return false;
  else {
    buff[i]--; // decrease
    i++;
    while (i < arrange_size) {
      buff[i] = total_size[i] - 1;   // reset the tail
      i++;
    }
    return true;
  }
}

Arrangement &Arrangement::operator=(Arrangement const &a){
  if(arrange_size) 
    delete buff;

  allocate(a.total_size, a.arrange_size);
  for(unsigned i = 0; i < arrange_size; i++)
    buff[i] = a.buff[i];
  return *this;  
}


bool Arrangement::operator==(Arrangement const &a){
  if((arrange_size != a.arrange_size) || 
     memcmp(total_size, a.total_size, arrange_size) != 0)
    return false;
  for(unsigned i = 0; i < arrange_size; i++)
    if(buff[i] != a.buff[i]) 
      return false;
  return true;  
}


bool Arrangement::operator<(Arrangement const &a){
  if((arrange_size != a.arrange_size) ||
     memcmp(total_size, a.total_size, arrange_size) != 0)
    std::cerr << "operator<() ... different dimension" << std::endl;
  for(unsigned i = 0; i < arrange_size; i++)
    if(buff[i] >= a.buff[i]) 
      return false;
  return true;  
}

