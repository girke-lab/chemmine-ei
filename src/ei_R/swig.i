%{
#include "desc.h"
#include "simpledb.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include "db_build.h"
#include <vector>
#include <iterator>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "solver.h"


int db2db_distance2file(char *dbFile, char *iddb1File,char *iddb2File,char *outfile);
int db2db_distance2file(char *dbFile, char *db2File,char *outfile);
int binaryCoord(char *inFile,char *outFile,int dim);

Solver* getSolver(int r,int d, double *refCoords);
//void embedCoord(Solver *s,int d,double *dist,double *result);
//double* embedCoord(Solver *s,int d,double *dist);
%}

%include "carrays.i"
%array_class(double,doubleArray);
int batch_sdf_parse(const char* sdfile, const char* dbfile);


int db2db_distance2file(char *dbFile, char *iddb1File,char *iddb2File,char *outfile);
int db2db_distance2file(char *dbFile, char *db2File,char *outfile);
int binaryCoord(char *inFile,char *outFile,int dim);


Solver* getSolver(int r,int d, double *refCoords);
//void embedCoord(Solver *s,int d,double *dist,double *result);
//double* embedCoord(Solver *s,int d,double *dist);
