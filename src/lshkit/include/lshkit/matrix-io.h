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

#ifndef LSHKIT_MATRIX_IO
#define LSHKIT_MATRIX_IO

#include <cassert>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>

namespace lshkit{

template <class T>
void Matrix<T>::load (const std::string &path)
{
  unsigned header[3]; /* entry size, row, col */
  assert(sizeof header == 3*4);
  int fd = ::open(path.c_str(), O_RDONLY);
  assert(fd >= 0);
  if (::read(fd, header, sizeof header) != sizeof header) assert(0);
  assert(header[0] == sizeof(float));
  reset(header[2], header[1]);
  unsigned long sz = sizeof(T) * dim * N;
	std::cerr << sz << " bytes to read into matrix." << std::endl;
  unsigned long r_sz;
  unsigned long r_total = 0;
  unsigned long file_sz = sz;
  do {
  r_sz = ::read(fd, (char*) dims + r_total, sz);
  sz -= r_sz;
  r_total += r_sz;
  } while (sz > 0 and r_sz != 0);
  if (sz != 0) {
    std::cerr << "matrix file error: read " << r_total << "bytes, when "
              << file_sz << " bytes are mentioned in the header." << std::endl;
    assert(0);
  }
  ::close(fd);
}

template <class T>
void Matrix<T>::dump (const std::string &path)
{
  unsigned header[3];
  int fd = ::open(path.c_str(), O_WRONLY | O_CREAT, 0644);
  assert(fd >= 0);
    header[0] = sizeof(float);
    header[1] = N;
    header[2] = dim;
    if (::write(fd, header, sizeof header) != sizeof header) assert(0);
  unsigned long sz = sizeof(T) * dim * N;
  unsigned long w_sz;
  unsigned long w_total = 0;
  unsigned long matrix_sz = sz;
  do {
  w_sz = ::write(fd, (char*) dims + w_total, sz);
  sz -= w_sz;
  w_total += w_sz;
  } while (sz > 0 and w_sz != 0);
  if (sz != 0) {
    std::cerr << "matrix file error: rite " << w_total << "bytes, when "
              << matrix_sz << " bytes are in the matrix." << std::endl;
    assert(0);
  }
  if (::write(fd, dims, sz) != sz) assert(0);
  ::close(fd);
}


}

#endif

