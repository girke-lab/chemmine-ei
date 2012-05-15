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

#include <fstream>

namespace lshkit {

/* matrix */
template <class T>
class Matrix
{
	unsigned int dim, N;
	T *dims;
	T **vecs;

	void load (const char *);
	void dump (const char *);
public:
	void reset (unsigned int _dim, unsigned int _N)
	{
		dim = _dim;
		N = _N;
		if (dims != NULL) delete[] dims;
		if (vecs != NULL) delete[] vecs;
		size_t _s = dim;
		_s *= N;
		dims = new T[_s];
		vecs = new T*[N];
		for (unsigned int i = 0; i < N; i++) {
			vecs[i] = dims + i * dim;
		}
	}
	void free (void) {
		dim = N = 0;
		if (dims != NULL) delete[] dims;
		if (vecs != NULL) delete[] dims;
		dims = NULL;
		vecs = NULL;
	}

	Matrix () :dim(0), N(0), dims(NULL), vecs(NULL) {}
	Matrix (unsigned int _dim, unsigned int _N) : dims(NULL), vecs(NULL) { reset(_dim, _N); }
	~Matrix () { if (dims != NULL) delete[] dims; if (vecs != NULL) delete[] vecs; }

	const T *operator [] (unsigned int i) const { return vecs[i]; }
	T *operator [] (unsigned int i) { return vecs[i]; }
	unsigned int getDim () const {return dim; }
	unsigned int getSize () const {return N; }

	void load (const std::string &path);
	void dump (const std::string &path);

	Matrix (std::string path): dims(NULL),vecs(NULL) { load(path); }
};

}

#include <lshkit/matrix-io.h>
