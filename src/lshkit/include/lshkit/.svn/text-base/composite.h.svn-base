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

#ifndef __LSHKIT_COMPOSITE__
#define __LSHKIT_COMPOSITE__

/*
   \file composite.h
   \brief Defining a bunch of LSH composite operations.
 
   LSH composition is to use an existing LSH class as building
   block to generate a new class of LSH.
*/

namespace lshkit {

//  The modulo operation on hash values.
/*
    The modulo of an LSH function by some value N is usually still locality
    sensitive.  This can be used to limit the hash value of certain LSH, so that
    the hash value can be used to index a fixed-sized hash table.

    Let LSH be the original class, and N be the divisor, then the parameter type
    is defined as

        struct Parameter {
            unsigned first;     // the divisor
            LSH::Parameter second;  // the parameter to the original LSH
        };

    The domain is the same as the original LSH.

    The range is always N.
*/
template <typename LSH>
class Tail
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));
protected:
    LSH lsh_;
    unsigned range_;
public:
    typedef std::pair<unsigned, typename LSH::Parameter> Parameter;
    typedef typename LSH::Domain Domain;

    Tail()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        range_ = param.first;
        lsh_.reset(param.second, rng);
    }

    template <typename RNG>
    Tail(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return range_;
    }

    unsigned operator () (Domain obj) const
    {
        return lsh_(obj) % range_;
    }
};

// The modulo operation with fixed divisor.

/*
   This is the same as Tail, except the divisor is determined at compile time
    and passed in as a template parameter.  Both the Domain and Parameter of the
    original LSH are kept the same.
*/

template <typename LSH, unsigned range>
class FixedTail
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));
protected:
    LSH lsh_;
public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;

    FixedTail()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        lsh_.reset(param, rng);
    }

    template <typename RNG>
    FixedTail(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return range;
    }

    unsigned operator () (Domain obj) const
    {
        return lsh_(obj) % range;
    }

    LSH &getLsh ()
    {
        return lsh_;
    }
};

// Take the Least Significant Bit of the hash value.
/*
   This is a special case of FixedTail, with the divisor being 2.
*/
template <typename LSH>
class LSB: public FixedTail<LSH, 2>
{
public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;

    LSB () {}

    template <typename RNG>
    LSB(const Parameter &param, RNG &rng)
        : FixedTail<LSH, 2>(param, rng)
    {
    }
};

// LSB for DeltaLSHes.
/*
   The original LSH must be a DeltaLSH.
 */
template <typename LSH>
class DeltaLSB: public FixedTail<LSH, 2>
{
    BOOST_CONCEPT_ASSERT((DeltaLshConcept<LSH>));

    typedef FixedTail<LSH, 2> Super;

public:
    typedef typename LSH::Parameter Parameter;
    typedef typename LSH::Domain Domain;

    DeltaLSB() {}

    template <typename RNG>
    DeltaLSB(const Parameter &param, RNG &rng)
        : Super(param, rng)
    {
    }
    unsigned operator () (Domain obj, float *delta) const
    {
        float d;
        unsigned r = Super::getLsh()(obj, &d);
        delta = min(d, 1.0F - d);
        return r;
    }
};

// Concatenation of a number of LSHes.
/*
   The concatenation of a number of LSHes of the same class is usually used as a
   new LSH to augment the locality sensitiveness.  LSH composition implements
   such functionality to concatenate N independent LSHes.

   The domain of the original LSH is kept the same.  The new parameter is is
   defined as:
        struct Parameter {
            unsigned first;     // # of LSHes to concatenate
            LSH::Parameter second;  // the parameter to the original LSH
        };

   Because the hash value is represented as unsigned int, which has only 32
   bits, the range of the original LSH need to be small enough so that the
   concatenated value does not overflow.  Specifically, we require that
                           LSH::getRange()^N <= 2^32.
   We depend on the user to guarantee that the range of LSH only depends on the
   parameter, so an array of such LSHes initialized with the same parameter but
   independent random numbers have the same range.
 */
template <typename LSH>
class Repeat
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));

protected:
    std::vector<LSH> lsh_;
    unsigned dup_;
    unsigned range_;
    unsigned unit_;

public:
    typedef std::pair<unsigned, typename LSH::Parameter> Parameter;
    typedef typename LSH::Domain Domain;

    Repeat ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        dup_ = param.first;
        assert(dup_ > 0);
        lsh_.resize(dup_);
        lsh_[0].reset(param.second, rng);
        range_ = unit_ = lsh_[0].getRange();
        assert(unit_ > 0);
        assert( (unsigned)(1 << (sizeof(unsigned) * 8 / dup_)) >= unit_);
        for (unsigned i = 1; i < dup_; ++i)
        {
            lsh_[i].reset(param.second, rng);
            assert(unit_ == lsh_[i].getRange());
            range_ *= unit_;
        }
    }

    template <typename RNG>
    Repeat(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return range_;
    }

    unsigned operator () (Domain obj) const
    {
        unsigned ret = 0;
        for (unsigned i = 0; i < dup_; ++i)
        {
            ret *= unit_;
            ret += lsh_[i](obj);
        }
        return ret;
    }
};


// Apply a random hash to the concatenation a number of hash values.
/*
   This composition is to workaround the case where the range of individual
   LSHes are so large that the concatenation can not be held in a single
   unsgined int.  The method is to further hash the concatenated value.
   Specifically, if <h1, h2, ..., hN> are the original values, this composition
   produces (a1*h1 + a2*h2 + aN*hN), with a1~aN being random unsigned integers
   sampled at random.  The range of the produced LSh is 0.

   The domain of the original LSH is kept the same.  The new parameter is is
   defined as:
        struct Parameter {
            unsigned first;     // # of LSHes to concatenate
            LSH::Parameter second;  // the parameter to the original LSH
        };
 */

template <typename LSH>
class RepeatHash
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));

protected:
    std::vector<LSH> lsh_;
    std::vector<unsigned> a_;

public:
    typedef std::pair<unsigned, typename LSH::Parameter> Parameter;
    typedef typename LSH::Domain Domain;

    RepeatHash ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        assert(param.first > 0);
        lsh_.resize(param.first);
        for (unsigned i = 0; i < param.first; ++i)
        {
            lsh_[i].reset(param.second, rng);
        }
        a_.resize(param.first);

        for (unsigned i = 0; i < param.first; ++i) a_[i] = rng();
    }

    template <typename RNG>
    RepeatHash(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return 0;
    }

    unsigned operator () (Domain obj) const
    {
        unsigned ret = 0;
        for (unsigned i = 0; i < lsh_.size(); ++i)
        {
            ret += lsh_[i](obj) * a_[i];
        }
        return ret;
    }
};


//  XOR of a number of 1-bit LSHes.
/*
   The XOR of a number of 1-bit LSHes has higher locality sensitiveness than
   the original LSHes.  This serves a similar purpose as RepeatHash.

   The domain of the original LSH is kept the same.  The new parameter is is
   defined as:
        struct Parameter {
            unsigned first;     // # of LSHes to XOR
            LSH::Parameter second;  // the parameter to the original LSH
        };
*/
template <typename LSH>
class Xor
{
    BOOST_CONCEPT_ASSERT((LshConcept<LSH>));
    std::vector<LSH> lsh_;
public:
    typedef std::pair<unsigned, typename LSH::Parameter> Parameter;
    typedef typename LSH::Domain Domain;

    Xor ()
    {
    }

    template <typename RNG>
    void reset(const Parameter &param, RNG &rng)
    {
        lsh_.resize(param.first);
        for (unsigned i = 0; i < lsh_.size; ++i)
        {
            lsh_[i].reset(param.second, rng);
            assert(lsh_[i].getRange() == 2);
        }
    }

    template <typename RNG>
    Xor(const Parameter &param, RNG &rng)
    {
        reset(param, rng);
    }

    unsigned getRange () const
    {
        return 2;
    }

    unsigned operator () (Domain obj) const
    {
        unsigned ret = 0;
        for (unsigned i = 0; i < lsh_.size(); ++i)
        {
            ret = ret ^ lsh_[i](obj);
        }
        return ret;
    }
};

// Delta version of XOR.
/*
   This is essentially the same as XOR.  The delta of XOR is the
   minimal of the deltas of all the original LSHes.
*/
template <typename LSH>
class DeltaXor: public Xor<LSH>
{
public:
    typedef typename Xor<LSH>::Parameter Parameter;
    typedef typename Xor<LSH>::Domain Domain;

    DeltaXor ()
    {
    }

    template <typename RNG>
    DeltaXor(const Parameter &param, RNG &rng): Xor<LSH>(param, rng)
    {
    }

    unsigned operator () (Domain obj, float *delta) const
    {
        unsigned ret = 0;
        float m = std::numeric_limits<float>::max(), d;
        for (unsigned i = 0; i < Xor<LSH>::lsh_.size(); ++i)
        {
            ret = ret ^ Xor<LSH>::lsh_[i](obj, &d);
            m = min(m, d);
        }
        *delta = m;
        return ret;
    }
};

}

#endif

