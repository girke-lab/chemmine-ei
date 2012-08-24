

#include <R.h>
#include <Rinternals.h>

#include "solver.h"

extern "C" {

   SEXP embedCoord(SEXP s, SEXP d, SEXP dist);
   SEXP embedCoordTest(SEXP r,SEXP d, SEXP refCoords, SEXP dist);
}

// refcoords is r x d
Solver* getSolver(int r,int d, double *refCoords)
{

   Solver *s = new Solver(d,r,3,refCoords);
   
   return s;
}
SEXP embedCoordTest(SEXP r,SEXP d, SEXP refCoords, SEXP dist)
{
   
   SEXP ans;

   int i,j;
   //for(i=0; i< INTEGER(r)[0]; i++)
    //  for(j=0;j<INTEGER(d)[0];j++)
     //    Rprintf("refCoord[%d][%d]=%f\n",i,j,REAL(refCoords)[i+INTEGER(r)[0]*j]);


   //Rprintf("before pointer cast\n");
   Solver s(INTEGER(d)[0],INTEGER(r)[0],3,REAL(refCoords));
   //Rprintf("before allocation\n");

   //for(i=0;i<INTEGER(d)[0];i++)
      //Rprintf("%f ",REAL(dist)[i]);
   //Rprintf("\n");
   
   PROTECT(ans = allocVector(REALSXP, INTEGER(d)[0]));
   //Rprintf("before optim\n");
   s.optim(REAL(ans),REAL(dist));
   //Rprintf("before unprotect\n");
   UNPROTECT(1);
   //Rprintf("done\n");
   return ans;
}


SEXP embedCoord(SEXP s, SEXP d, SEXP dist)
{
   
   SEXP ans;
   //Rprintf("before pointer cast\n");
   Solver *solver = reinterpret_cast<Solver*>( R_ExternalPtrAddr(s));
   //Rprintf("before allocation\n");
   PROTECT(ans = allocVector(REALSXP, INTEGER(d)[0]));
   //Rprintf("before optim\n");
   solver->optim(REAL(ans),REAL(dist));
   //Rprintf("before unprotect\n");
   UNPROTECT(1);
   //Rprintf("done\n");
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
