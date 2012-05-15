
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/dynamic_bitset.hpp>
#include <lshkit/mplsh.h>
#include <lshkit/matrix.h>
#include <lshkit/eval.h>
#include <sys/time.h> 
#include <vector>
#include <fstream>
#include <queue>
#include <iostream>

#ifndef _NO_REFINE_
#include "../../atompair/desc.h"
#include "../../atompair/simpledb.h"
#endif

class Timer
{
	struct	timeval	start; 
public:
    Timer () {}
    void tick ()
    {
	    gettimeofday(&start, 0); 
    }
    void tuck (const char *msg) const
    {
        struct timeval end;
	    float	diff; 
	    gettimeofday(&end, 0); 

	    diff = (end.tv_sec - start.tv_sec) 
	   			+ (end.tv_usec - start.tv_usec) * 0.000001; 
        std::cerr << msg << ':' <<  diff << std::endl;
    }
};


using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 

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

int main (int argc, char *argv[])
{
    string data_file, query_file, cdb_file;
//    string benchmark;

    float W;
    unsigned M, L, H;
    unsigned K, T;
    unsigned nseg, seg;
	unsigned int k2 = 0;
	bool do_refine = false;

    Timer timer;
	std::vector<unsigned int> queries;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message.")
		(",W", po::value<float>(&W)->default_value(1.282), "")
		(",M", po::value<unsigned>(&M)->default_value(11), "")
		(",L", po::value<unsigned>(&L)->default_value(10), "")
		(",H", po::value<unsigned>(&H)->default_value(1017881), "")
		(",K", po::value<unsigned>(&K)->default_value(10), "")
		(",T", po::value<unsigned>(&T)->default_value(10), "")
		("data,D", po::value<string>(&data_file), "")
		("shrink,S", po::value<unsigned int>(&k2)->default_value(0), "")
		("cdb,C", po::value<string>(&cdb_file), "")
//		("benchmark,B", po::value<string>(&benchmark), "")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);	

	if (vm.count("help") || (vm.count("data") < 1))
	{
		cerr << desc;
		return 0;
	}

#ifndef _NO_REFINE_
	//IndexedDB db;
	PreloadedDB db;
	std::fstream ifs;
	if (vm.count("cdb") >= 1) {
		cerr << "using cdb to refine" << endl;
		if (not db.open(cdb_file.c_str())) return 1;
		cerr << "cdb loaded" << endl;
		do_refine = true;
	}
	if (do_refine and (k2 == 0 or k2 > K)) {
		cerr << "set k2 to 1/5 of K" << endl;
		k2 = K / 5;
	}
#endif

    cerr << "Loading data...";
    Matrix<float> data(data_file);
    cerr << "done." << endl;

    /*
    Benchmark<> bench(K, Q);
    cout << "Loading benchmark...";
    bench.load(benchmark);
    cout << "done." << endl;

    for (unsigned i = 0; i < Q; ++i)
    {
        for (unsigned j = 0; j < K; ++j)
        {
            assert(bench.getAnswer(i)[j].key < data.getSize());
				}
    }
    */

    cerr << "Initializing index..." << endl;

    typedef MultiProbeLshIndex<MatrixAccessor> Index;
        
    Index::Parameter param;

    param.W = W;
    param.H = H;
    param.M = M;
    param.dim = data.getDim();
    DefaultRng rng;

    MatrixAccessor accessor(data);

    Index index(param, rng, accessor, L);

    cerr << "done." << endl;

    cerr << "Populating index..." << endl;

    timer.tick();

    {
        boost::progress_display progress(data.getSize(), cerr);
        for (unsigned i = 0; i < data.getSize(); ++i)
        {
            index.insert(i);
            ++progress;
        }
    }
    timer.tuck("CREATE");

    cerr << "Running queries..." << endl;


	std::cerr << "Type query index number. type 0 to exit" << std::endl;
	while (true)
	{
		int q_id;
		std::cerr << ">>";
		std::cin >> q_id;
		if (! std::cin.good()) {
			if (std::cin.eof())
				break;
			std::cin.clear();
			std::cin.ignore(1024, '\n');
			std::cerr << "Invalid input. Try it again." << std::endl;
			continue;
		}

		if (q_id == 0) break;
		q_id --; 

		if (q_id >= data.getSize()) {
			std::cerr << "Input is too large. Maximum is "
				      << data.getSize() - 1 
					  <<" Try it again." << std::endl;
			continue;
		}

		timer.tick();
		unsigned cnt;
		Topk<unsigned> topk;
		topk.reset(K);

		topk.reset(K);
		index.query(data[q_id], topk, T, &cnt);
#ifndef _NO_REFINE_		
		if (do_refine) {
			std::vector<unsigned int> desc_query;
			db.at(q_id, desc_query);
			std::vector<unsigned int> dbcmp;
			std::priority_queue<std::pair<double, unsigned int> > q_dist;
			for (unsigned int j = 0; j < K; j ++) {
				unsigned int hit = topk[j].key;
				db.at(hit, dbcmp);
				double s = similarity(desc_query, dbcmp);
				if (q_dist.size() == k2) {
					if (1-s > q_dist.top().first) continue;
					q_dist.pop();
					q_dist.push(std::pair<double, unsigned int>(1-s, hit));
				} else {
					q_dist.push(std::pair<double, unsigned int>(1-s, hit));
				}
			}
			std::vector<std::pair<double, unsigned int> > r(k2);
			int idx = 0;
			while (not q_dist.empty()) {
				r[idx ++] = q_dist.top();
				q_dist.pop();
			}
			for (int i = 0; i < k2; i ++)
				std::cout << r[k2-i-1].second + 1<< ":" << r[k2-i-1].first << " ";
		} else {
#endif

		for (unsigned j = 0; j < K; j ++)
			cout << topk[j].key + 1 << ":" << topk[j].dist << " ";

#ifndef _NO_REFINE_		
		}
#endif
		cout << endl;
		timer.tuck("Query time");
	}

    return 0;
}


// vim:tabstop=4:shiftwidth=4:smartindent
