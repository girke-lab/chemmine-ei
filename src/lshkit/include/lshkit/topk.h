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


#ifndef __LSHKIT_TOPK__
#define __LSHKIT_TOPK__

#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>

namespace lshkit {

template <typename KEY>
struct TopkEntry
{
    KEY key;
    float dist;   
    bool match (const TopkEntry &e) const { return key == e.key; }
    bool match (KEY e) const { return key == e; }

    TopkEntry (KEY key_, float dist_) : key(key_), dist(dist_) {}
    TopkEntry () { dist = std::numeric_limits<float>::max(); }
    void reset () { dist = std::numeric_limits<float>::max(); }

    friend bool operator < (const TopkEntry &e1, const TopkEntry &e2)
    {
        return e1.dist < e2.dist;
    }
};

template <class KEY>
class Topk: public std::vector<TopkEntry<KEY> >
{
public:
    typedef TopkEntry<KEY> Element;
    typedef typename std::vector<TopkEntry<KEY> > Base;

    Topk () {}

	Topk (unsigned k): Base(k) {}

	~Topk () {}

	void reset (unsigned k) { this->resize(k); for (typename Base::iterator it = this->begin(); it != this->end(); ++it) it->reset(); }

    const Element &back () const { return Base::back(); }

    Topk &operator << (Element t)
	{
        unsigned i = this->size() - 1;
        unsigned j;
        if (t < this->at(i))
        {
            for (;;)
            {
                if (i == 0) break;
                j = i - 1;
                if (this->at(j).match(t)) return *this; 
                if (this->at(j) < t) break;
                i = j;
            }
            /* i is the place to insert to */

            j = this->size() - 1;
            for (;;)
            {
                if (j == i) break;
                this->at(j) = this->at(j-1);
                --j;
            }
            this->at(i) = t;
        }
        return *this;
	}

    float recall (const Topk<KEY> &topk /* to be evaluated */) const
    {
        unsigned matched = 0;
        for (typename Base::const_iterator ii = this->begin(); ii != this->end(); ++ii)
        {
            for (typename Base::const_iterator jj = topk.begin(); jj != topk.end(); ++jj)
            {
                if (ii->match(*jj))
                {
                    matched++;
                    break;
                }
            }
        }
        return float(matched)/float(this->size());
    }
};


};

#endif

