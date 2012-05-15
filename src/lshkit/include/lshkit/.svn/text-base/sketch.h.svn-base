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

#ifndef __LSHKIT_SKETCH__
#define __LSHKIT_SKETCH__

/*
   LSH-based Sketch Construction

   Sketches are compact representations (as bit-vectors) of large objects.
   The distance between the original objects can be approximated by the hamming
   distance between sketches.

   For more information on sketches and asymmetric distance estimators, see

       Wei Dong, Moses Charikar, Kai Li. Asymmetric Distance Estimation with
       Sketches for Similarity Search in High-Dimensional Spaces . To appear in
       Proceedings of the 31st Annual International ACM SIGIR Conference on Research
       & Development on Information Retrieval. Singapore. July 2008.
*/

namespace lshkit {

/*
   This class uses 1-bit LSH to generate sketches.  The sketches are stored as
   arrays of type CHUNK (unsigned char by default).
*/
template <typename LSH, typename CHUNK = unsigned char>
class Sketch
{
    BOOST_CONCEPT_ASSERT((DeltaLshConcept<LSH>));

    unsigned dim_;
    std::vector<LSH> lsh_;

public:
    typedef typename LSH::Parameter Parameter; // LSH parameter
    typedef typename LSH::Domain Domain; // Domain of LSH & Sketch
    static const unsigned CHUNK_BIT = sizeof(CHUNK) * 8; // #bits in CHUNK

    /* Constructor
        bits    - #bits in the sketch
        param   - parameter to LSH
        engine  - random number generator
     */
    template <typename Engine>
    Sketch(unsigned bits, Parameter param, Engine &engine) : 
        dim_((bits + CHUNK_BIT - 1) / CHUNK_BIT), lsh_(dim_ * CHUNK_BIT)
    {
        for (unsigned i = 0; i < lsh_.size(); ++i) lsh_[i].reset(param, engine);
    }

    /* sketch construction
       out shoud point to an array of ceil(nbit/CHUNK_BIT) elements.
       */
	void apply (Domain in, CHUNK *out) const
	{
        unsigned l = 0;
		for (unsigned i = 0; i < dim_; ++i) {
			out[i] = 0;
			for (unsigned j = 0; j < CHUNK_BIT; ++j) {
                unsigned k = lsh_[l++](in);
				out[i] = out[i] | (k << j);
			}
		}
	}
    
    /* sketch construction
       also calculate the values used for asymmetric distance estimation.
     */
    void apply (Domain in, CHUNK *out, float *asym) const
    {
        unsigned l = 0;
		for (unsigned i = 0; i < dim_; ++i) {
			out[i] = 0;
			for (unsigned j = 0; j < CHUNK_BIT; ++j) {
                unsigned k = lsh_[l](in, &asym[l]);
				out[i] = out[i] | (k << j);
                l++;
			}
		}
    }

};


}

#endif

