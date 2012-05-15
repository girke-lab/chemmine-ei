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


#ifndef __LSHKIT_FLAT__
#define __LSHKIT_FLAT__

#include <algorithm>
#include <lshkit/topk.h>

namespace lshkit {

template <typename LSH, typename ACCESSOR, typename METRIC>
class FlatIndex
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));
public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;
    typedef typename ACCESSOR::Key Key;


protected:
    typedef std::vector<Key> Bin;
    std::vector<LSH> lshs_;
    std::vector<std::vector<Bin> > tables_;
    ACCESSOR accessor_;
    METRIC metric_;

public:
    template <typename Engine>
    FlatIndex(const Parameter &param, Engine &engine, const ACCESSOR &accessor,
            const METRIC metric, unsigned L) : 
        lshs_(L), tables_(L), accessor_(accessor), metric_(metric)
    {
        for (unsigned i = 0; i < L; ++i)
        {
            lshs_[i].reset(param, engine);
            tables_[i].resize(lshs_[i].getRange());
        }
    }

    void insert (Key key)
    {
        for (unsigned i = 0; i < lshs_.size(); ++i)
        {
            unsigned index = lshs_[i](accessor_(key));
            tables_[i][index].push_back(key);
        }
    }

    void query (const Domain &obj, Topk<Key> &topk, unsigned *pcnt = (unsigned *)0)
    {
//        if (L == 0) L = lshs_.size();
 //       assert(L <= lshs_.size());
        unsigned L = lshs_.size();
        unsigned cnt = 0;
        for (unsigned i = 0; i < L; ++i)
        {
            unsigned index = lshs_[i](obj);
            Bin &bin = tables_[i][index];
            for (typename Bin::const_iterator it = bin.begin();
                    it != bin.end(); ++it)
            {
                ++cnt;
                topk << typename Topk<Key>::Element(*it, metric_(obj, accessor_(*it)));
            }
        }
        if (*pcnt != 0) *pcnt = cnt;
    }
};


}

#endif

