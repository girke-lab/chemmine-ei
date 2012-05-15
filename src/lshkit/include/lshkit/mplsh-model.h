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
#ifndef __LSHKIT_MPLSH_MODEL__
#define __LSHKIT_MPLSH_MODEL__

#include <fstream>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace lshkit {

typedef boost::math::normal_distribution<double> NormalDist;
typedef boost::math::gamma_distribution<double> GammaDist;

GammaDist GammaDistMLE (double M, double G);

class DataParam
{
    double M, G;
    double a_M, b_M, c_M;
    double a_G, b_G, c_G;
    /*
    DataParam () {}
    DataParam (const DataParam &dp)
        : M(dp.M), G(dp.G),
        a_M(dp.a_M), b_M(dp.b_M), c_M(dp.c_M),
        a_G(dp.a_G), b_G(dp.b_G), c_G(dp.c_G) {}
        */
public:
    DataParam (const std::string &path)
    {
        std::ifstream is(path.c_str());
        is >> M >> G >> a_M >> b_M >> c_M
             >> a_G >> b_G >> c_G;
    }

    GammaDist globalDist () const
    {
        return GammaDistMLE(M, G);
    }

    GammaDist topkDist (unsigned N, unsigned K) const
    {
        double m, g;
        m = std::exp(a_M) * std::pow(double(N), b_M) * std::pow(double(K), c_M);
        g = std::exp(a_G) * std::pow(double(N), b_G) * std::pow(double(K), c_G);
        return GammaDistMLE(m,g);
    }
};

class MultiProbeLshModel
{
    unsigned L_;
    double W_;
    unsigned M_, T_;
public:
    MultiProbeLshModel(unsigned L, double W, unsigned M, unsigned T)
        : L_(L), W_(W), M_(M), T_(T)
    {}

    double recall (double l2dist) const;

    void setL (unsigned L) { L_ = L; }
    void setW (double W) { W_ = W; }
    void setT (unsigned T) { T_ = T; }
    void setM (unsigned M) { M_ = M; }

    unsigned getT () const {return T_; }

};

class MultiProbeLshDataModel: public MultiProbeLshModel
{
// temporary, not thread safe
    GammaDist globalDist_;
    std::vector<GammaDist> topkDists_;

public:
    MultiProbeLshDataModel(const DataParam &param, unsigned N, unsigned K)
        : MultiProbeLshModel(0,0,0,0), globalDist_(1.0),
        topkDists_(K, globalDist_)
    {
        globalDist_ = param.globalDist();
        for (unsigned k = 0; k < K; k++)
        {
            topkDists_[k] = param.topkDist(N, k + 1);
        }
    }

    double avgRecall () const;
    double cost () const;
};

class MultiProbeLshRecallTable
{
	unsigned step_;
	float min_, max_;
    boost::numeric::ublas::matrix<float> table_;
public:

    MultiProbeLshRecallTable (MultiProbeLshModel model, unsigned d_step, float d_min,
            float d_max)
        : step_(d_step), min_(d_min), max_(d_max), table_(model.getT(), d_step)
    {
        unsigned T = model.getT();
        double delta = (max_ - min_) / step_;
        for (unsigned t = 0; t < T; ++t)
        {
            model.setT(t+1);
            for (unsigned d = 0; d < step_; ++d)
            {
                table_(t, d) = model.recall(sqr(min_ + delta * d));
            }
        }
    }
    
    float lookup (float dist, int T)
    {
        unsigned d;
        if (dist < min_) return 1.0;
        if (!(dist < max_)) return 0.0;
        d = (unsigned int) ((dist - min_) * step_ / (max_ - min_));
        return table_(T-1,d);
    }

};

}

#endif
