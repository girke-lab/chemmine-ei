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
#include <lshkit/common.h>
#include <lshkit/mplsh.h>
#include <lshkit/mplsh-model.h>
#include <btune.h>

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 

#define MIN_W	0.4
#define MAX_W	20
#define NUM_W	10000
#define DELTA_W	((MAX_W - MIN_W) / NUM_W)


// T, L, M, W
struct btune_minmax minmax[]= {{1, 100},{1, 20},{4, 20},{0, NUM_W}};

double target_recall;

int fun (double *f, const int *x, void *_model)
{
    MultiProbeLshDataModel *model = reinterpret_cast<MultiProbeLshDataModel *>(_model);
	model->setW(MIN_W + DELTA_W * x[3]);
	model->setM(x[2]);
	model->setL(x[1]);
    model->setT(x[0]);

	f[0] = model->cost();
	f[1] = model->avgRecall() - target_recall;
	return 0;
}

int x[] = {0, -1, -1, -1};
double y[2];

int main (int argc, char *argv[])
{
	double w;
    unsigned N, K;
    string data_param;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message.")
		(",T", po::value<int>(&x[0])->default_value(-1), "")
		(",L", po::value<int>(&x[1])->default_value(-1), "")
		(",M", po::value<int>(&x[2])->default_value(-1), "")
		(",W", po::value<double>(&w)->default_value(-1), "")
		("num,N", po::value<unsigned>(&N), "")
		("data,D", po::value<string>(&data_param), "")
		("recall,R", po::value<double>(&target_recall)->default_value(0.9), "")
		("topk,K", po::value<unsigned>(&K)->default_value(20), "")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);	

	if (vm.count("help") || (vm.count("data") < 1) || (vm.count("recall") < 1) ||
            (vm.count("topk") < 1) || (vm.count("num") < 1))
	{
		cout << desc;
		return 0;
	}

	if (w <= 0) x[3] = -1;
	else x[3] = (w - MIN_W) / DELTA_W;

    MultiProbeLshDataModel model(DataParam(data_param), N, K);

	struct btune btune;
	btune_init(&btune, 4, minmax, fun, &model);
	btune_tune(&btune, x, y);
	btune_cleanup(&btune);

	printf("T = %d\tL = %d\tM = %d\tW = %g\trecall = %g\tcost = %g\n",
		x[0], x[1], x[2], MIN_W + DELTA_W * x[3], target_recall + y[1], y[0]);

	return 0;
}

