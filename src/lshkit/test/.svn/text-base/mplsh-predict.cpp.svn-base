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

using namespace std;
using namespace lshkit;
namespace po = boost::program_options; 

int main (int argc, char *argv[])
{
	double w;
    unsigned N, K, T, L, M;
    string data_param;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message.")
		(",T", po::value<unsigned>(&T)->default_value(1), "")
		(",L", po::value<unsigned>(&L)->default_value(1), "")
		(",M", po::value<unsigned>(&M)->default_value(1), "")
		(",W", po::value<double>(&w)->default_value(1), "")
		("num,N", po::value<unsigned>(&N), "")
		("data,D", po::value<string>(&data_param), "")
		("topk,K", po::value<unsigned>(&K)->default_value(20), "")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);	

	if (vm.count("help") || (vm.count("data") < 1) || 
            (vm.count("topk") < 1) || (vm.count("num") < 1))
	{
		cout << desc;
		return 0;
	}

    MultiProbeLshDataModel model(DataParam(data_param), N, K);
    model.setT(T);
    model.setL(L);
    model.setM(M);
    model.setW(w);

    cout << model.avgRecall() << '\t' << model.cost() << endl;

	return 0;
}

