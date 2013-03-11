%{
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "solver.h"


int binaryCoord(char *inFile,char *outFile,int dim);
Solver* getSolver(int r,int d, double *refCoords);
int eucsearch2file(const char* matrix,const char* queryMatrix,int n_results,char* outfile);
%}

int binaryCoord(char *inFile,char *outFile,int dim);
Solver* getSolver(int r,int d, double *refCoords);
int eucsearch2file(const char* matrix,const char* queryMatrix,int n_results,char* outfile);
