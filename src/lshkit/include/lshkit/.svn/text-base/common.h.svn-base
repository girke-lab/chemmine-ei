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

#ifndef __LSHKIT_COMMON__
#define __LSHKIT_COMMON__


#include <cmath>
#include <limits>
#include <vector>
#include <boost/random.hpp>

// Concept checking requires boost > 1.35.
// remove the following line to disable concept checking
#define CONCEPT_CHECK 1

#ifdef CONCEPT_CHECK
#include <lshkit/concept.h>
#else
#define BOOST_CONCEPT_ASSERT(A)
#endif


namespace lshkit {

// The default random number generator.
typedef boost::mt19937 DefaultRng;

// Some of the frequently used distributions.
typedef boost::normal_distribution<float> Gaussian;
typedef boost::cauchy_distribution<float> Cauchy;
typedef boost::uniform_real<float> UniformReal;
typedef boost::uniform_int<int> UniformInt;
typedef boost::uniform_int<unsigned> UniformUnsigned;

/*
   The minimum & maximum of two values.
   We define this so that <algorithm> doesn't have to be included solely for
   this simple function.
 */
template <typename T>
T min (T a, T b) { return a < b ? a : b; }

template <typename T>
T max (T a, T b) { return a < b ? b : a; }

// Squre function.
template <typename T> T sqr (const T &x) { return x * x; }

}


#endif

