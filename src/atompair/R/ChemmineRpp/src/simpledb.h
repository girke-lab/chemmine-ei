#ifndef SIMPLEDB_H_
#define SIMPLEDB_H_
#include "debug.h"
#include <vector>
#include <fstream>
int serialize(std::vector<unsigned int> & descs, std::fstream & ofs);
int load(std::vector<unsigned int> & descs, std::fstream & ifs);
int load(unsigned int ** p_buf, std::fstream & ifs);

#define HEADER_SIZE 16

#endif /*SIMPLEDB_H_*/
