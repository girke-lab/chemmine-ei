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

#ifndef __LSHKIT_PROBE__
#define __LSHKIT_PROBE__

#include <vector>
#include <fstream>
#include <string>
#include <stdint.h>

#include <lshkit/common.h>
#include <lshkit/lsh.h>
#include <lshkit/composite.h>
#include <lshkit/metric.h>
#include <lshkit/flat-index.h>
#include <lshkit/mplsh-model.h>

namespace lshkit
{

// the data structure

struct Probe
{
    uint32_t mask;
    uint32_t shift;
    float score;
    unsigned reserve;
    bool operator < (const Probe &p) const { return score < p.score; }
    Probe operator + (const Probe &m) const
    {
        Probe ret;
        ret.mask = mask | m.mask;
        ret.shift = shift | m.shift;
        ret.score = score + m.score;
        return ret;
    }
    bool conflict (const Probe &m)
    {
        return (mask & m.mask) != 0;
    }
    static const unsigned MAX_M = 20;
    static const unsigned MAX_T = 100;
}; 

typedef std::vector<Probe> ProbeSequence;

void GenProbeSequenceTemplate (ProbeSequence &seq, unsigned M, unsigned T);

class ProbeSequenceTemplates: public std::vector<ProbeSequence>
{
public:
    ProbeSequenceTemplates(unsigned max_M, unsigned max_T)
        : std::vector<ProbeSequence>(max_M + 1)
    {
        for (unsigned i = 1; i <= max_M; ++i)
        {
            GenProbeSequenceTemplate(at(i), i, max_T);
        }
    }
};

extern ProbeSequenceTemplates __probeSequenceTemplates;


class MultiProbeLsh: public RepeatHash<GaussianLsh> 
{
    typedef RepeatHash<GaussianLsh> Super;
public:
    struct Parameter {
        unsigned H;
        unsigned M;
        unsigned dim;
        float W;
    };

    typedef Super::Domain Domain;
private:
    unsigned H_;


public:
    MultiProbeLsh () {}

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        Super::Parameter bp;
        bp.first = param.M;
        bp.second.dim = param.dim;
        bp.second.W = param.W;
        H_ = param.H;
        Super::reset(bp, rng);
    }

    template <typename RNG>
    MultiProbeLsh(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return H_;
    }

    unsigned operator () (Domain obj) const
    {
        return Super::operator ()(obj) % H_;
    }

    void genProbeSequence (Domain obj, std::vector<unsigned> &seq, unsigned T);

};


template <typename ACCESSOR>
class MultiProbeLshIndex: public FlatIndex<MultiProbeLsh, ACCESSOR, metric::l2<float> >
{
    typedef FlatIndex<MultiProbeLsh, ACCESSOR, metric::l2<float> > Super;

    MultiProbeLshModel model_;
    MultiProbeLshRecallTable recall_;
public: 
    typedef typename Super::Parameter Parameter;
    typedef typename Super::Domain Domain;
    typedef typename Super::Key Key;

    template <typename Engine>
    MultiProbeLshIndex(const Parameter &param, Engine &engine, const ACCESSOR &accessor, unsigned L) 
        : Super(param, engine, accessor, metric::l2<float>(param.dim), L),
        model_(L, param.W, param.M, Probe::MAX_T),
        recall_(model_, 200, 0.00001, 1.0)
    {
    }

    void query (const Domain &obj, Topk<Key> &topk, unsigned T, unsigned *pcnt = (unsigned *)0)
    {
        std::vector<unsigned> seq;
        unsigned L = Super::lshs_.size();
        unsigned cnt = 0;
        Super::accessor_.reset();
        for (unsigned i = 0; i < L; ++i)
        {
            Super::lshs_[i].genProbeSequence(obj, seq, T);
            for (unsigned j = 0; j < seq.size(); ++j)
            {
                typename Super::Bin &bin = Super::tables_[i][seq[j]];
                for (typename Super::Bin::const_iterator it = bin.begin();
                    it != bin.end(); ++it)
                    if (Super::accessor_.mark(*it))
                    {
                        ++cnt;
                        topk << typename Topk<Key>::Element(*it, Super::metric_(obj,
                                    Super::accessor_(*it)));
                    }
            }
        }
        if (pcnt != 0) *pcnt = cnt;
    }

    void query (const Domain &obj, Topk<Key> &topk, float recall, unsigned *pcnt = (unsigned *)0)
    {
        unsigned L = Super::lshs_.size();
        std::vector<std::vector<unsigned> > seqs(L);
        for (unsigned i = 0; i < L; ++i) Super::lshs_[i].genProbeSequence(obj,
                seqs[i], Probe::MAX_T);

        unsigned cnt = 0;
        Super::accessor_.reset();
        for (unsigned j = 0; j < Probe::MAX_T; ++j)
        {
            if (j >= seqs[0].size()) break;
            for (unsigned i = 0; i < L; ++i)
            {
                typename Super::Bin &bin = Super::tables_[i][seqs[i][j]];
                for (typename Super::Bin::const_iterator it = bin.begin();
                    it != bin.end(); ++it)
                    if (Super::accessor_.mark(*it))
                    {
                        ++cnt;
                        topk << typename Topk<Key>::Element(*it, Super::metric_(obj,
                                    Super::accessor_(*it)));
                    }
            }
            float r = 0.0;
            for (unsigned i = 0; i < topk.size(); ++i)
            {
                r += recall_.lookup(topk[i].dist, j+1);
            }
            r /= topk.size();
            if (r >= recall) break;
        }
        if (pcnt != 0) *pcnt = cnt;
    }
};

// model

}


#endif

