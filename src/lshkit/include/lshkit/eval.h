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

#ifndef __LSHKIT_EVAL__
#define __LSHKIT_EVAL__

#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>

namespace lshkit {

template <typename KEY = unsigned>
class Benchmark
{
	int K_, Q_;

    std::vector<unsigned> queries_;
    std::vector<Topk<KEY> > topks_;
public:
	Benchmark(): K_(0), Q_(0) {}
	/* do random initialization if max_index > 0*/
	void reset(int K, int Q) {
		K_ = K;
		Q_ = Q;
        queries_.resize(Q);
        topks_.resize(Q);
    	for (int i = 0; i < Q; i++) {
				topks_[i].reset(K);
		}
	}

	Benchmark(int K, int Q)
			 {reset(K, Q);}

	~Benchmark() {}

	void load (std::istream &is)
	{
		std::string line;
        unsigned k;
		for (int i = 0; i < Q_; i++) {
            is >> queries_[i];
            is >> k;
			for (int j = 0; j < K_; j++) {
				is >> topks_[i][j].key;
				is >> topks_[i][j].dist;
			}
			getline(is, line);
		}
	}
	void dump (std::ostream &os) const
	{
		for (int i = 0; i < Q_; i++) {
            os << queries_[i] << '\t' << K_ << '\t';
			for (int j = 0; j < K_; j++) {
				os << '\t' <<  topks_[i][j].key;
				os << '\t' <<  topks_[i][j].dist;
			}
			os << std::endl;
		}
	}

	void load (std::string path)
	{
		std::ifstream is(path.c_str());
		load(is);
		is.close();
	}
	void dump (std::string path) const
	{
		std::ofstream os(path.c_str());
		dump(os);
		os.close();
	}

    unsigned getK () const { return K_; }
    unsigned getQ () const { return Q_; }

    unsigned getQuery (unsigned n) const {  return queries_[n]; }
	const Topk<KEY> &getAnswer (unsigned n) const { return topks_[n]; }
};

class Stat
{
	int count;
	float sum;
	float sum2;
	float min;
	float max;
public:
	Stat () : count(0), sum(0), sum2(0), min(std::numeric_limits<float>::max()), max(-std::numeric_limits<float>::max()) {}
	~Stat () {}
	void append (float r)
	{
		count++;
		sum += r;
		sum2 += r*r;
		if (r > max) max = r;
		if (r < min) min = r;
	}
	Stat & operator<< (float r) { append(r); return *this; }
	int getCount() { return count; }
	float getSum() { return sum; }
	float getAvg() { return sum/count; }
	float getMax() { return max; }
	float getMin() { return min; }
	float getStd()
	{
		if (count > 1) return sqrt((sum2 - (sum/count) * sum)/(count - 1)); 
		else return 0; 
	}
	void merge (Stat &stat)
	{
		count += stat.count;
		sum += stat.sum;
		sum2 += stat.sum2;
		if (stat.min < min) min = stat.min;
		if (stat.max > max) max = stat.max;
	}
};

template <typename RNG>
void SampleQueries (std::vector<unsigned> *qry, unsigned Q, unsigned max, RNG
        rng)
{
    boost::variate_generator<RNG &, UniformUnsigned> gen(rng, UniformUnsigned(0,
                max-1));
    qry->resize(Q);
    for (unsigned i = 0; i < Q; ++i)
    {
		for (;;)
        {
            qry->at(i) = gen();
            unsigned j;
            for (j = 0; j < i; j++) if (qry->at(i) == qry->at(j)) break;
            if (j >= i) break;
        }
    }
}

}

#endif

