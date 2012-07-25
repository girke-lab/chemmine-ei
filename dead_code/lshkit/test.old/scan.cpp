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

#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <lshkit/common.h>
#include <lshkit/mplsh-file.h>
#include <lshkit/matrix.h>
#include <lshkit/metric.h>
#include <lshkit/eval.h>
#include <sys/time.h> 

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
        std::cout << msg << ':' <<  diff << std::endl;
    }
};

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 


int main (int argc, char *argv[])
{
    Timer timer;
    string data_file;
    string query_file;

    unsigned K, Q;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message.")
		(",Q", po::value<unsigned>(&Q)->default_value(1), "")
		(",K", po::value<unsigned>(&K)->default_value(1), "")
		("data,D", po::value<string>(&data_file), "")
		("benchmark,B", po::value<string>(&query_file), "")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);	

	if (vm.count("help") || (vm.count("data") < 1) || (vm.count("benchmark") < 1))
	{
		cout << desc;
		return 0;
	}

    timer.tick();
    Matrix<float> data(data_file);

    Benchmark<unsigned> bench(K, Q);
    bench.load(query_file);

    Topk<unsigned> topk;

    metric::l2<float> l2(data.getDim());

    timer.tick();
    for (unsigned i = 0; i < Q; ++i)
    {
        topk.reset(K);
        unsigned q = bench.getQuery(i);
        for (unsigned j = 0; j < data.getSize(); ++j)
        {
            topk << Topk<unsigned>::Element(j, l2(data[q],
                                    data[i]));
        }
    }

    timer.tuck("QUERY");

    return 0;
}

