

#include <R.h>
#include <Rinternals.h>

#include "solver.h"

extern "C" {

   SEXP embedCoord(SEXP s, SEXP d, SEXP dist);
}

// refcoords is r x d
Solver* getSolver(int r,int d, double *refCoords)
{

   Solver *s = new Solver(d,r,3,refCoords);
   
   return s;
}

SEXP embedCoord(SEXP s, SEXP d, SEXP dist)
{
   
   SEXP ans;
   Solver *solver = reinterpret_cast<Solver*>( R_ExternalPtrAddr(s));
   PROTECT(ans = allocVector(REALSXP, INTEGER(d)[0]));
   solver->optim(REAL(ans),REAL(dist));
   UNPROTECT(1);
   return ans;
}

/*
//dist is length r, result is length d
//void embedCoord(Solver *s,int d,double *dist,double *result)
double* embedCoord(Solver *s,int d,double *dist)
{
   //s->optim(result,dist);
   double *embedded = new double[d];
   s->optim(embedded,dist);
   return embedded;
}
*/
