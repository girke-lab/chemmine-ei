/* 
    Copyright (C) 2008 Wei Dong <wdong@princeton.edu>. All Rights Reserved.
  
    This file is part of LSHKIT.
  
    LSHKIT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LSHKIT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LSHKIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <gsl/gsl_multifit.h>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <lshkit/common.h>
#include <lshkit/matrix.h>
#include <lshkit/topk.h>
#include <lshkit/eval.h>
#include <lshkit/metric.h>

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 

int main (int argc, char *argv[])
{
    string data_file;
    unsigned N, Q, K, F;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message.")
		(",N", po::value<unsigned>(&N)->default_value(50000), "")
		(",Q", po::value<unsigned>(&Q)->default_value(1000), "")
		(",K", po::value<unsigned>(&K)->default_value(100), "")
		(",F", po::value<unsigned>(&F)->default_value(10), "")
		("data,D", po::value<string>(&data_file), "")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);	

	if (vm.count("help") || (vm.count("data") < 1))
	{
		cout << desc;
		return 0;
	}

	/* load matrix */
    Matrix<float> data(data_file);

    metric::l2sqr<float> l2(data.getDim());

    DefaultRng rng;
    boost::variate_generator<DefaultRng &, UniformUnsigned> gen(rng,
             UniformUnsigned(0, data.getSize()-1));


	double gM = 0.0;
    double gG = 0.0;
    {
	    for (unsigned k = 0; k < N; ++k)
        {
            double dist, logdist;
            for (;;)
            {
                unsigned i = gen() % data.getSize();
                unsigned j = gen() % data.getSize();
                dist = l2(data[i], data[j]);
                logdist = log(dist);
                if (isnormal(logdist)) break;
            }
            gM += dist;
            gG += logdist;
        }
    	gM /= N;
    	gG /= N;
	    gG = exp(gG);
    }
    cout << gM << '\t' << gG << endl;

	/* sample query */
    vector<unsigned> qry;


    SampleQueries(&qry, Q, data.getSize(), rng);

	/* do the queries */
    vector<Topk<unsigned> > topks(Q);
    for (unsigned i = 0; i < Q; ++i) topks[i].reset(K);

/* ... */
	gsl_matrix *X = gsl_matrix_alloc(F * K, 3);
	gsl_vector *yM = gsl_vector_alloc(F * K);
	gsl_vector *yG = gsl_vector_alloc(F * K);
	gsl_vector *pM = gsl_vector_alloc(3);
	gsl_vector *pG = gsl_vector_alloc(3);
	gsl_matrix *cov = gsl_matrix_alloc(3,3);

    vector<double> M(K);
    vector<double> G(K);

    boost::progress_display progress(F, cerr);
    unsigned m = 0;
	for (unsigned l = 0; l < F; l++)
	{
        // Scan
        for (unsigned i = l; i< data.getSize(); i += F)
        {
            for (unsigned j = 0; j < Q; j++) 
            {
                int id = qry[j];
                if (i != id) 
                {
                    float d = l2(data[id], data[i]);
                    if (isnormal(d))
                    topks[j] << Topk<unsigned>::Element(i, d);
                }
            }
        }

        fill(M.begin(), M.end(), 0.0);
        fill(G.begin(), G.end(), 0.0);

		for (unsigned i = 0; i < Q; i++)
		{
			for (unsigned k = 0; k < K; k++)
			{
				M[k] += topks[i][k].dist;
				G[k] += log(topks[i][k].dist);
			}
		}

		for (unsigned k = 0; k < K; k++)
		{
			M[k] = log(M[k]/Q);
			G[k] /= Q;
			gsl_matrix_set(X, m, 0, 1.0);
			gsl_matrix_set(X, m, 1, log(double(data.getSize() * (l + 1)) / double(F)));
			gsl_matrix_set(X, m, 2, log(double(k + 1)));
			gsl_vector_set(yM, m, M[k]);
			gsl_vector_set(yG, m, G[k]);
            ++m;
		}

        ++progress;
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(F * K, 3);

    double chisq;

	gsl_multifit_linear(X, yM, pM, cov, &chisq, work);
	gsl_multifit_linear(X, yG, pG, cov, &chisq, work);

	cout << gsl_vector_get(pM, 0) << '\t'
		 << gsl_vector_get(pM, 1) << '\t'
	     <<	gsl_vector_get(pM, 2) << endl;
	cout << gsl_vector_get(pG, 0) << '\t'
		 << gsl_vector_get(pG, 1) << '\t'
	     <<	gsl_vector_get(pG, 2) << endl;

	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(pM);
	gsl_vector_free(pG);
	gsl_vector_free(yM);
	gsl_vector_free(yG);

	return 0;
}

