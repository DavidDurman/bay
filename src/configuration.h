/**
 * configuration.h
 * @description ADT for all possible value assignments for a node's configuration template, i.e. configuration.
 * @author David Durman, 2008
 */

#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

#include <cstring>
#include <iostream>

class configuration {
public:
  unsigned* buff;	// value assignments
  unsigned size;	// number of parents + me
  unsigned node;	// node id

  configuration(unsigned* b, unsigned s, unsigned n){
    node = n;
    size = s;
    buff = new unsigned[s];
    for (unsigned i = 0; i < s; i++)
      buff[i] = b[i];
    //    buff = (unsigned*) memcpy(buff, b, s * sizeof(unsigned));
  }

  bool operator==(const configuration& c) const {
    if (size != c.size || node != c.node)
      return false;
    return (bool)(memcmp(buff, c.buff, c.size * sizeof(unsigned)) == 0);
  }

  /**
   * Convert to decimal number (using Horner schema) and compare.
   */
  bool operator<(const configuration& c) const {
    if (node < c.node)
      return true;
    if (node > c.node)
      return false;
    if (size < c.size)
      return true;
    if (size > c.size)
      return false;

    int i1 = 0, i2 = 0;
    for (unsigned i = 0; i < size; i++){
      i1 = buff[i] + i1 * 10;
      i2 = c.buff[i] + i2 * 10;
    }
    if (i1 < i2)
      return true;
    return false;
  }

  void clear(){ delete [] buff; }

  void print(std::ostream& os){
    if (size == 0) return;
    os << node << ":";	// node
    os << buff[0];	// first
    for (unsigned i = 1; i < size; i++)
      os << " " << buff[i];
  }
};

#endif // _CONFIGURATION_H
