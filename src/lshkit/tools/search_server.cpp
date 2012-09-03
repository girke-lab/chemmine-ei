// derived from isearch.cpp
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

#include "../../atompair/profiling.h"
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/format.hpp>
#include <lshkit.h>
#include <vector>
#include <fstream>
#include <queue>
#include <string.h>
#include <stdlib.h>
#include <signal.h>

/**
  * \file isearch.cpp
  * \brief Interactive search using LSH index
  *
  * This program is an example of using MPLSH index.
  *
  * The program reconstruct the LSH index by default.  You can give the
  * --index option to make the program save the LSH index.  The next
  * time you run the program with the same --index option, the program
  * will try to load the previously saved index.  When a saved index is
  * used, you need to make sure that the dataset and other parameters match
  * the previous run.  However, the benchmark file, Q and K can be different.
  *
\verbatim
Allowed options:
  -h [ --help ]                   produce help message.
  -W [ -- ] arg (=1)
  -M [ -- ] arg (=1)
  -T [ -- ] arg (=1)              # probes
  -L [ -- ] arg (=1)              # hash tables
  -K [ -- ] arg (=50)             # nearest neighbor to retrieve
  -R [ -- ] arg (=3.40282347e+38) R-NN distance range
  -D [ --data ] arg               data file
  --index arg                     index file
  -H [ -- ] arg (=1017881)        hash table size, use the default value.
\endverbatim
  */

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 
const unsigned int BUFSIZE = 10485760;
const char* INPUT_FP = "coord.in";

/*
    You must provide an access class to query the MPLSH.
    MPLSH only saves keys (pointers to the real feature vectors) in the
    hash tables and it will rely on the accessor class to retrieve
    the feature vector.

    An accessor must provide three methods:

        bool mark (unsigned key);

        mark that key has been accessed.  If key has already been marked, return false,
        otherwise return true.  MPLSH will use this to avoid scanning the data more than
        one time per query.

        void reset ();

        to clear all the marks.
        
        const float *operator () (unsigned key);

        given a key, return the pointer to a feature vector.

    The MatrixAccessor class is used to access feature vectors stored in a Matrix.
*/

/* This class has been merged to include/lshkit/matrix.h */
/*
class MatrixAccessor
{
    const Matrix<float> &matrix_;
    boost::dynamic_bitset<> flags_;
public:
    typedef unsigned Key;
    MatrixAccessor(const Matrix<float> &matrix)
        : matrix_(matrix), flags_(matrix.getSize()) {}

    MatrixAccessor(const MatrixAccessor &accessor)
        : matrix_(accessor.matrix_), flags_(matrix_.getSize()) {}

    void reset ()
    {
        flags_.reset();
    }

    bool mark (unsigned key)
    {
        if (flags_[key]) return false;
        flags_.set(key);
        return true;
    }

    const float *operator () (unsigned key)
    {
        return matrix_[key];
    }
};
*/

void on_sig(int status) 
{
	std::cerr << "waked up to work: " << status << std::endl;
	return; 
}

int main (int argc, char *argv[])
{

		signal(SIGINT, on_sig);
    string data_file;
    string index_file;

		Timer t;
    float W, R, desired_recall = 1.0;
    unsigned M, L, H;
    unsigned Q, K, T;
    bool use_index = false; // load the index from a file
    bool do_refine = false; // whether to perform refinement step
    unsigned skip = 0;

    vector<unsigned> queries;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message.")
        (",W", po::value<float>(&W)->default_value(1.0), "")
        (",M", po::value<unsigned>(&M)->default_value(1), "")
        (",T", po::value<unsigned>(&T)->default_value(1), "# probes")
        (",L", po::value<unsigned>(&L)->default_value(1), "# hash tables")
        (",K", po::value<unsigned>(&K)->default_value(0), "default # nearest neighbor to retrieve")
        ("radius,R", po::value<float>(&R)->default_value(numeric_limits<float>::max()), "R-NN distance range (L2)")
        ("data,D", po::value<string>(&data_file), "data file")
        ("index", po::value<string>(&index_file), "index file")
        (",H", po::value<unsigned>(&H)->default_value(1017881), "hash table size, use the default value.")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm); 

    if (vm.count("help") || (vm.count("data") < 1))
    {
        cerr << desc;
        return 0;
    }

    if (vm.count("radius") >= 1) {
        R *= R; // we use L2sqr in the program.
    }

    if (vm.count("index") == 1) {
        use_index = true;
    }

    cerr << "LOADING DATA..." << endl;
    FloatMatrix data(data_file);

    typedef MultiProbeLshIndex<unsigned> Index;

    FloatMatrix::Accessor accessor(data);
    Index index;

    // try loading index
    bool index_loaded = false;

    // load index?
    if (use_index) {
        ifstream is(index_file.c_str(), ios_base::binary);
        if (is) {
            is.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);
            cerr << "LOADING INDEX..." << endl;
            index.load(is);
            //verify(is);
            index_loaded = true;
        }
    }

    if (!index_loaded) {
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
        cerr << "CONSTRUCTING INDEX..." << endl;

        {
            boost::progress_display progress(data.getSize(), cerr);
            for (unsigned i = 0; i < data.getSize(); ++i)
            {
                // Insert an item to the hash table.
                // Note that only the key is passed in here.
                // MPLSH will get the feature from the accessor.
                index.insert(i, data[i]);
                ++progress;
            }
        }

        if (use_index) {
            cerr << "SAVING INDEX..." << endl;
            {
                ofstream os(index_file.c_str(), ios_base::binary);
                os.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);
                index.save(os);
            }
        }
    }

    metric::l2sqr<float> l2sqr(data.getDim());
		char* inp = new char[BUFSIZE];
		float* query_vec = new float[data.getDim()];
		while (true) {
			unsigned int k = K;
			std::cout << ">>";
			std::cout.flush();
			pause();
			ifstream inp_f;
			inp_f.open(INPUT_FP, std::ios::in);
			if (not inp_f.good()) {
				std::cerr << "Error in opening input file\n" << std::endl;
				continue;
			}
			inp_f.getline(inp, BUFSIZE);
			if (inp_f.fail()) {
				std::cerr << "Error in reading input file\n" << std::endl;
				continue;
			}

			// parsing input line
			char* str = strtok(inp, " ");
			if (strcmp(str, "/Q") == 0) break;
			if (strcmp(str, "/K") == 0) {
					str = strtok(NULL, " ");
					int _k = atoi(str);
					if (_k > 0) k = _k;
					str = strtok(NULL, " ");
			}
			bool vector_ok = true;
			for (unsigned dim = 0; dim < data.getDim(); dim ++) {
				query_vec[dim] = atof(str);
				if (dim < data.getDim() - 1) {
					str = strtok(NULL, " ");
					if (str == NULL) {
						std::cerr << "Wrong number of dimensions in input vector" 
							<< std::endl;
						vector_ok = false;
						break;
					}
				}
			}
			if (not vector_ok) continue;

         unsigned cnt;
         Topk<unsigned> topk;
         float maxValue = std::numeric_limits<float>::max();
         TopkScanner<FloatMatrix::Accessor, metric::l2sqr<float> > query(accessor, l2sqr, k, R);
			topk.reset(k);
			query.reset(query_vec);
			t.start();

			//cout<<"topk init:"<<endl;
			//for (unsigned j = 0; j < k; j ++)
			//	cout << topk[j].key + 1 << ":" << topk[j].dist << " ";
			//cout<<endl;

			//cout<<"query vector: "<<endl;

			//for (unsigned j = 0; j < data.getDim(); j ++)
			//	cout << query_vec[j] << " ";
			//cout<<endl;


			index.query(query_vec, T, query);
			topk.swap(query.topk());
  
			cout << "OK:"; 
			for (unsigned j = 0; j < k; j ++)
				if(topk[j].dist != maxValue)
					cout << topk[j].key + 1 << ":" << topk[j].dist << " ";
  
			cout << endl;
			cerr << boost::format("QUERY TIME: %1%s.") % t.pause() << endl;
			t.reset();
    }

		if (inp) {
			delete[] inp;
			inp = NULL;
		}
    return 0;
}

