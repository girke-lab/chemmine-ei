#include <lshkit.h>
#include <vector>
#include <queue>
#include <R.h>
#include <Rinternals.h>

using namespace lshkit;

extern "C" {
   SEXP lshsearch(SEXP queries, SEXP matrixFile, SEXP Rin );
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
// add: W,M,T,L,K,H
SEXP lshsearch(SEXP queries, SEXP matrixFile, SEXP Rin )
{
   float W=1.0;
   unsigned H = 1017881;
   unsigned M = 1;
   unsigned L = 1;
   unsigned K = 600;
   unsigned T = 1;

   //float R = REAL(Rin)[0] >= 1? REAL(Rin)[0] * REAL(Rin)[0] : REAL(Rin)[0];
   float R = std::numeric_limits<float>::max();

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
   float *queryPtr=(float *)REAL(queries);
   double *queryPtr2=REAL(queries);
   

   //for(int i=0; i < numQueries; i++)
   for(int i=0; i < 1; i++)
   {
      Rprintf("query %d:\n",i);
      for(int j=0;j<querySize;j++)
         //Rprintf("%f ",REAL(queries)[j]);
         Rprintf("%f ",queryPtr[j]);
      Rprintf("\n");
      continue;

      unsigned cnt;
      Topk<unsigned> topk;
      float maxValue = std::numeric_limits<float>::max();
      TopkScanner<FloatMatrix::Accessor, metric::l2sqr<float> > query(accessor, l2sqr, K, R);
      topk.reset(K);

      query.reset(queryPtr);
      
      index.query(queryPtr, T, query);
      topk.swap(query.topk());
      queryPtr+=querySize;


      for (unsigned j = 0; j < K; j ++)
         Rprintf("%d:%f ",topk[j].key,topk[j].dist);
      Rprintf("\n");

      for(unsigned j = 0; j < K; j++)
      {
         REAL(result)[j*2+i*numQueries] = topk[j].key;
         if(topk[j].dist != maxValue)
            REAL(result)[j*2+i*numQueries+1] = topk[j].dist;
         else
            REAL(result)[j*2+i*numQueries+1] = 0;
      }

   }
   UNPROTECT(1);
   return result;
}

