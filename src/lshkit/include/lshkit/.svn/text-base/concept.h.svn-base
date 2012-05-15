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

#ifndef __LSHKIT_CONCEPT__
#define __LSHKIT_CONCEPT__

#include <boost/concept/assert.hpp>
#include <boost/concept/usage.hpp>

namespace lshkit {

/*
    This class checks the LSH concept.  A LSH class should define the following
    items to be used with LSHKIT:

    typedef ... Parameter;  // The parameter type.
    typedef ... Domain;     // The domain of the LSH function.
    LSH();                  // Default constructor.
    void reset(const Parameter &, RNG &);   // Initialize the object
    LSH(const Parameter &, RNG &);  // equivalent to LSH() followed by reset()
    unsigned getRange () const; // get the range of the hash value
    unsigned operator () const (const Domain);  // apply the hash function.

    The type Parameter should contain only POD, so that it can be written to
    /load from files in the internal binary representation.

    If getRange() returns 0, then the hash value (returned by operator()) can be
    any thing.  Otherwise, it can be from 0 to getRange() - 1.

    The reset() function and on of the constructors take a random number
    generator.  This random number generator should not be further used after
    reset() or the constructor has returned.
*/

template <typename LSH>
struct LshConcept
{
    typedef typename LSH::Domain Domain;
    typedef typename LSH::Parameter Parameter;

    BOOST_CONCEPT_USAGE(LshConcept)
    {
        LSH lsh1(param, rng1);
        LSH lsh2(param, rng2);
        LSH lsh3;
        lsh.reset(param, rng1);
        lsh.reset(param, rng2);
        same_type(lsh.getRange(), u);
        same_type(lsh(object), u);
    }
protected:
    boost::mt11213b rng1;
    boost::mt19937 rng2;
    LSH lsh;
    Domain object;
    const Parameter param;

    unsigned u;
    template <typename T>
    void same_type(T const &, T const&);
};

/*
   Some of the LSH functions are created by rouding a real number to an integer,
   and the rounded off part (delta) of the hash value carries information that
   sometimes can be useful.  Such LSH functions are represented by DeltaLSH,
   which should implement the following besides all that apply to the general
    LSH:

    unsigned operator () const (Domain, float *delta);  // apply the hash function.
*/
template <typename LSH>
struct DeltaLshConcept: public LshConcept<LSH>
{
    typedef LshConcept<LSH> Super;
    BOOST_CONCEPT_USAGE(DeltaLshConcept)
    {
        float delta;
        Super::same_type(Super::lsh(Super::object, &delta), Super::u);
    }
};

}

#endif

