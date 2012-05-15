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


#ifndef __LSHKIT_HISTOGRAM__
#define __LSHKIT_HISTOGRAM__

namespace lshkit {

template <typename LSH>
class Histogram
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));

    std::vector<LSH> lsh_;

    unsigned M_, N_;
    unsigned dim_;
    unsigned unit_;
public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;

    template <typename RNG>
    Histogram(unsigned M, unsigned N, Parameter parameter, RNG &rng)
        : lsh_(M * N), M_(M), N_(N)
    {
        for (unsigned i = 0; i < lsh_.size(); ++i)
        {
            lsh_[i].reset(parameter, rng);
        }
        unit_ = lsh_[0].getRange();
        dim_ = N_ * unit_;
    }

	int getDim () const
    {
        return dim_;
    }

    void zero (float *out)
    {
        std::fill(out, out + dim_, 0);
    }
	
	void add (float *out, Domain in, float weight = 1.0)
    {
        unsigned k = 0;
        unsigned base = 0;
        for (unsigned i = 0; i < N_; ++i)
        {
            for (unsigned j = 0; j < M_; ++j)
            {
                unsigned index = base + (lsh_[k])(in);
                out[index] +=  weight;
                ++k;
            }
            base += unit_;
        }
    }
};

}

#endif
