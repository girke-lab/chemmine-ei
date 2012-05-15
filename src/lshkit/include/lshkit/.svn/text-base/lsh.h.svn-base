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

#ifndef __LSHKIT_LSH__
#define __LSHKIT_LSH__

// This file defines a set of basic LSH families.

namespace lshkit {

/*
   Stable distribution based LSH

   This LSH is defined on the D-dimensional vector space.  For a vector X, the
   hash value is defined as
                  h(X) = [b + a1*X1 + a2*X2 + ... + aD*XD] / W
   where W is a positive value called the window size; b is sampled uniformly in
   [0, W); a1 ~ aD are N random variables sampled from the so-called stable
   distribution, specifiled by the template parameter DIST,  independently at
   random.

   The domain of the LSH is (const float *), and the parameter is defined as
        struct Parameter {
            unsigned dim;
            float W;
        };
   The range of this LSH is 0.

   Following are two spacial cases:
       typedef StableDistLsh<Cauchy> CauchyLsh;
       typedef StableDistLsh<Gaussian> GaussianLsh;
   Cauchy distribution is 1-stable and Gaussian distribution is 2-stable.  These
   two LSHes can be used to approximate L1 and L2 distances respectively.
   
   For more information on stable distribution based LSH, see

       Mayur Datar , Nicole Immorlica , Piotr Indyk , Vahab S. Mirrokni,
       Locality-sensitive hashing scheme based on p-stable distributions,
       Proceedings of the twentieth annual symposium on Computational geometry, June
       08-11, 2004, Brooklyn, New York, USA.

*/
template <typename DIST>
class StableDistLsh
{
    std::vector<float> a_;
    float b_;
    float W_;
    unsigned dim_;
public:
    struct Parameter
    {
        unsigned dim;
        float W;
    };

    typedef const float *Domain;

    StableDistLsh ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        a_.resize(param.dim);
        W_ = param.W;
        dim_ = param.dim;

        boost::variate_generator<RNG &, DIST> gen(rng, DIST());

        for (unsigned i = 0; i < dim_; ++i) a_[i] = gen();

        b_ = boost::variate_generator<RNG &, UniformReal>(rng, UniformReal(0,W_))();
    }

    template <typename RNG>
    StableDistLsh(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }


    unsigned getRange () const
    {
        return 0;
    }

    unsigned operator () (Domain obj) const
    {
        float ret = b_;
        for (unsigned i = 0; i < dim_; ++i)
        {
            ret += a_[i] * obj[i];
        }
        return unsigned(int(std::floor(ret / W_)));
    }

    unsigned operator () (Domain obj, float *delta) const
    {
        float ret = b_;
        for (unsigned i = 0; i < dim_; ++i)
        {
            ret += a_[i] * obj[i];
        }
        ret /= W_;

        float flr =  std::floor(ret);
        *delta = ret - flr;
        return unsigned(int(flr));
    }
};

typedef StableDistLsh<Cauchy> CauchyLsh;
typedef StableDistLsh<Gaussian> GaussianLsh;

class HyperPlaneLsh
{
    std::vector<float> a_;

public:
    struct Parameter
    {
        unsigned dim;
    };
    typedef const float *Domain;

    HyperPlaneLsh ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        a_ = boost::variate_generator<RNG &, boost::uniform_on_sphere<float> >
                (rng, boost::uniform_on_sphere<float>(param))();
    }

    template <typename RNG>
    HyperPlaneLsh(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }


    unsigned getRange () const
    {
        return 2;
    }

    unsigned operator () (Domain obj) const
    {
        float ret = 0;
        for (unsigned i = 0; i < a_.size(); ++i)
        {
            ret += a_[i] * obj[i];
        }
        return ret >= 0 ? 1 : 0;
    }

    unsigned operator () (Domain obj, float *delta) const
    {
        float ret = 0;
        for (unsigned i = 0; i < a_.size(); ++i)
        {
            ret += a_[i] * obj[i];
        }
        if (ret >= 0)
        {
            *delta = ret;
            return 1;
        }
        else
        {
            *delta = -ret;
            return 0;
        }
    }
};

class ThresholdingLsh
{
    unsigned dim_;
    float threshold_;
public:
    struct Parameter
    {
        unsigned dim;
        float min;
        float max;
    };
    typedef const float *Domain;

    ThresholdingLsh ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        dim_ = boost::variate_generator<RNG &, UniformUnsigned>(rng, UniformUnsigned(0,param.dim - 1))();
        threshold_ = boost::variate_generator<RNG &, UniformReal>(rng, UniformReal(param.min,param.max))();
    }

    template <typename RNG>
    ThresholdingLsh(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }


    unsigned getRange () const
    {
        return 2;
    }

    unsigned operator () (Domain obj) const
    {
        return obj[dim_] >= 0 ? 1 : 0;
    }

    unsigned operator () (Domain obj, float *delta) const
    {
        float ret = obj[dim_] - threshold_;
        if (ret >= 0)
        {
            *delta = ret;
            return 1;
        }
        else
        {
            *delta = -ret;
            return 0;
        }
    }
};

}

#endif

