#include <lshkit.h>
#include <vector>
#include <queue>
#include <R.h>
#include <Rinternals.h>

using namespace lshkit;

extern "C" {
   SEXP lshsearch(SEXP queries, SEXP matrixFile, 
      SEXP Win, 
      SEXP Hin, 
      SEXP Min, 
      SEXP Lin, 
      SEXP Kin, 
      SEXP Tin, 
      SEXP Rin
   );
}

typedef MultiProbeLshIndex<unsigned> Index;

int loadIndex(Index &index, FloatMatrix &data, 
      float W, unsigned H, unsigned M, unsigned L)
{
   // We define a short name for the MPLSH index.
   Index::Parameter param;

   // Setup the parameters.  Note that L is not provided here.
   param.W = W;
   param.range = H; // See H in the program parameters.  You can just use the default value.
   param.repeat = M;
   param.dim = data.getDim();
   DefaultRng rng;

   index.init(param, rng, L);
   // The accessor.

   // Initialize the index structure.  Note L is passed here.


   for (unsigned i = 0; i < data.getSize(); ++i)
   {
       // Insert an item to the hash table.
       // Note that only the key is passed in here.
       // MPLSH will get the feature from the accessor.
       index.insert(i, data[i]);
   }


   /*
   if (use_index) {
      cerr << "SAVING INDEX..." << endl;
      {
          ofstream os(index_file.c_str(), ios_base::binary);
          os.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);
          index.save(os);
      }
   }
   */
}
int check(SEXP in,int deflt) {
   return INTEGER(in)[0] == NA_INTEGER ? deflt : INTEGER(in)[0];
}
double check(SEXP in,double deflt) {
   return ISNA(REAL(in)[0]) ? deflt : REAL(in)[0];
}
SEXP lshsearch(SEXP queries, SEXP matrixFile, 
      SEXP Win, 
      SEXP Hin, 
      SEXP Min, 
      SEXP Lin, 
      SEXP Kin, 
      SEXP Tin, 
      SEXP Rin
   )
{
   float W = check(Win,1.0);
   unsigned H = check(Hin, 1017881 );
   unsigned M = check(Min,1);
   unsigned L = check(Lin,1);
   unsigned K = check(Kin,600);
   unsigned T = check(Tin,1);
   float R = ISNA(REAL(Rin)[0])? std::numeric_limits<float>::max() : 
                                 (float)(REAL(Rin)[0]*REAL(Rin)[0]);
   Rprintf("W: %f H:%d M:%d L:%d K:%d T:%d R:%f\n",W,H,M,L,K,T,R);

   FloatMatrix data(CHAR(STRING_ELT(matrixFile,0)));


   FloatMatrix::Accessor accessor(data);
   Index index;

   loadIndex(index,data,W,H,M,L);

   metric::l2sqr<float> l2sqr(data.getDim());

   SEXP queryDim = getAttrib(queries,R_DimSymbol);
   int numQueries = INTEGER(queryDim)[1];
   int querySize = INTEGER(queryDim)[0];
   Rprintf("numQueries: %d, querySize: %d\n",numQueries,querySize);
   SEXP result;
   PROTECT(result = alloc3DArray(REALSXP,numQueries,K,2));
   float *queryPtr = new float[querySize];
   

   int k=0;
   for(int i=0; i < numQueries; i++)
   {
      //Rprintf("query %d:\n",i);
      for(int j=0;j<querySize;j++){
         queryPtr[j]=(float)REAL(queries)[k++];
         //Rprintf("%f ",queryPtr[j]);
      }
      //Rprintf("\n");

      unsigned cnt;
      Topk<unsigned> topk;
      float maxValue = std::numeric_limits<float>::max();
      TopkScanner<FloatMatrix::Accessor, metric::l2sqr<float> > query(accessor, l2sqr, K, R);
      topk.reset(K);

      query.reset(queryPtr);
      
      index.query(queryPtr, T, query);
      topk.swap(query.topk());


      //for (unsigned j = 0; j < K; j ++)
      //   if(topk[j].dist != maxValue)
      //      Rprintf("%d:%f ",topk[j].key,topk[j].dist);
      //Rprintf("\n");

      for(unsigned j = 0; j < K; j++)
      {
         int index=j*numQueries+i;  //for key
         int index2=(K+j)*numQueries+i; // for value
         if(topk[j].dist != maxValue){
            REAL(result)[index]   = topk[j].key+1;
            REAL(result)[index2] = topk[j].dist;
         }else{
            REAL(result)[index]   = -1;
            REAL(result)[index2] = -1;
         }
      }

   }
   delete queryPtr;
   UNPROTECT(1);
   return result;
}

