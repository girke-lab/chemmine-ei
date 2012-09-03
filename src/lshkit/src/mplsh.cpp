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

#include <cassert>
#include <queue>
#include <lshkit/mplsh.h>

#include <stdlib.h>

namespace lshkit
{
    ProbeSequenceTemplates __probeSequenceTemplates(Probe::MAX_M, Probe::MAX_T);

    void GenExpectScores (ProbeSequence &seq, unsigned M)
    {
        assert(M <= sizeof(seq[0].mask)* 8);
        seq.resize(2 * M);
        for (unsigned l = 0; l < M; ++l)
        {
            unsigned r = 2 * M - l - 1;
            seq[l].mask = seq[r].mask = seq[r].shift = leftshift(l);
            seq[l].shift = 0;
            seq[l].reserve = seq[r].reserve = 0;
            float delta = (l + 1.0) / (M + 1.0) * 0.5;
            seq[l].score = (l + 1.0) * (l + 2.0) / (M + 1.0) / (M + 2.0) * 0.25;
            seq[r].score = 1.0 - 2.0 * delta + seq[l].score;
        }
    }

    class ProbeGT
    {
    public:
        bool operator () (const Probe &p1, const Probe &p2) const
        {
            return p2 < p1;
        }
    };

    void GenProbeSequenceTemplate (ProbeSequence &seq, unsigned M, unsigned T)
    {
		printf("mark 1\n");
        ProbeSequence scores;
        GenExpectScores(scores, M);
        assert(T > 0);

		printf("mark 2\n");
        std::priority_queue<Probe, std::vector<Probe>, ProbeGT> heap;
        Probe init;
        init.mask = init.shift = 0;
        init.score = 0;
        init.reserve = 0;
		printf("mark 3\n");
        heap.push(init);
		printf("mark 4\n");

        seq.clear();
		printf("mark 5\n");

        for (;;)
        {
		printf("mark 51\n");
            if (heap.empty()) break;
		printf("mark 52\n");
            seq.push_back(heap.top());
		printf("mark 53\n");
            if (seq.size() == T) break;
		printf("mark 54\n");

            Probe shift = heap.top();
		printf("mark 541\n");
            heap.pop();
		printf("mark 55\n");

            for (unsigned next = shift.reserve; next < 2 * M; ++next)
            {
		printf("mark 56\n");
                if (!shift.conflict(scores[next]))
                {
                    Probe tmp = shift + scores[next];
                    tmp.reserve = next + 1;
                    heap.push(tmp);
                }
            }
        }
		printf("mark 6\n");
    }

    void MultiProbeLsh::genProbeSequence (Domain obj, std::vector<unsigned>
            &seq, unsigned T) const
    {
        ProbeSequence scores;
        std::vector<unsigned> base;
        scores.resize(2 * lsh_.size());
        base.resize(lsh_.size());
        for (unsigned i = 0; i < lsh_.size(); ++i)
        {
            float delta;
            base[i] = Super::lsh_[i](obj, &delta);
            scores[2*i].mask = i;
            scores[2*i].reserve = 1;    // direction
            scores[2*i].score = delta;
            scores[2*i+1].mask = i;
            scores[2*i+1].reserve = unsigned(-1);
            scores[2*i+1].score = 1.0 - delta;
        }
        std::sort(scores.begin(), scores.end());

        ProbeSequence &tmpl = __probeSequenceTemplates[lsh_.size()];

        seq.clear();
        for (ProbeSequence::const_iterator it = tmpl.begin();
                it != tmpl.end(); ++it)
        {
            if (seq.size() == T) break;
            const Probe &probe = *it;
            unsigned hash = 0;
            for (unsigned i = 0; i < lsh_.size(); ++i)
            {
                unsigned h = base[scores[i].mask];
                if (probe.mask & leftshift(i))
                {
                    if (probe.shift & leftshift(i))
                    {
                        h += scores[i].reserve;
                    }
                    else
                    {
                        h += unsigned(-1) * scores[i].reserve;
                    }
                }
                hash += h * a_[scores[i].mask];
            }
            seq.push_back(hash % H_);
        }
    }
}

