
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

#ifndef _NO_REFINE_
#include "desc.h"
#include "simpledb.h"
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
    string data_file, query_file, cdb_file, o_resume;
//    string benchmark;

    float W;
    unsigned M, L, H;
    unsigned K, T;
    unsigned nseg, seg;
	unsigned int k2 = 0;
	bool do_refine = false;
	unsigned int skip = 0;

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
		("query,Q", po::value<string>(&query_file), "")
		("shrink,S", po::value<unsigned int>(&k2)->default_value(0), "")
		("cdb,C", po::value<string>(&cdb_file), "")
		("nsegments,", po::value<unsigned>(&nseg), "")
		("segment,", po::value<unsigned>(&seg), "")
		("resume,", po::value<string>(&o_resume), "")
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
	IndexedDB db;
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

	if (vm.count("resume") >= 1) {
		cerr << "continue from " << o_resume << endl;
		std::ifstream fo_resume(o_resume.c_str());
		if (not fo_resume.good()) {
			cerr << "cannot read " << o_resume << endl;
			return 1;
		}
		char line[1048576];
		while (true) {
			fo_resume.getline(line, 1048576);
			if (fo_resume.good())
				skip ++;
			else
				break;
		}
		
		cerr << o_resume << " has " << skip << " lines." << endl;
	}

    cerr << "Loading data...";
    Matrix<float> data(data_file);
    cerr << "done." << endl;
    cerr << "Loading queries...";
	unsigned int n_queries;
	bool query_all = false;
	int offset = 0;
	if (vm.count("query") >= 1) {
		fstream qf(query_file.c_str(), std::ios::in);
		unsigned int temp;
		assert(qf.good());
		while (true) {
			qf >> temp;
			if (not qf.good()) break;
			queries.push_back(temp);
		}
		n_queries = queries.size();
		cerr << n_queries << " queries read" << endl;
	} else {
		if (vm.count("nsegments") < 1 || (vm.count("segment") < 1)) {
			cerr << desc;
			return 0;
		}
		int per_process = data.getSize() / nseg + 1;
		offset = per_process * seg;
		n_queries = per_process;
		if (n_queries + offset > data.getSize())
			n_queries = data.getSize() - offset;
		cerr << "Query [" << offset << " , " << n_queries + offset << ")" << endl;
		query_all = true;
	}

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


    timer.tick();
	boost::progress_display progress(n_queries, cerr);
	for (unsigned i = offset; i < offset + n_queries; ++i)
	{
		unsigned cnt;
		Topk<unsigned> topk;
		topk.reset(K);
		int q_id;
		if (not query_all) 
			q_id = queries[i];
		else
			q_id = i;

		if (skip > 0) {
			skip --;
		} else {

			while (true) {
				topk.reset(K);
				index.query(data[q_id], topk, T, &cnt);
	/*			bool pass = false;
				for (unsigned int j = 0; j < K; j ++) {
					unsigned int hit = topk[j].key;
					if (hit == q_id) {
						pass = true;
						break;
					}
				}
				if (pass) break;
				*/
				break;
			}
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
		}
		++progress;
	}
    timer.tuck("QUERY");

    return 0;
}


// vim:tabstop=4:shiftwidth=4:smartindent
