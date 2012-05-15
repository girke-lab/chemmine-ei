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

#ifndef __WDONG_BTUNE__
#define __WDONG_BTUNE__

#include <db.h>

struct btune_minmax
{
	int min;
	int max;	/* max will always work */
};

typedef int (*btune_f) (double *r, const int *p1, void *p2);

struct btune
{
	DB *db;
	int N;
	struct btune_minmax *minmax;
	btune_f fun;
	void *param;
};

/* min r[0]
   st
   	r[1] >= 0
*/

int btune_init (struct btune *btune, int N, const struct btune_minmax *minmax, btune_f f1, void *param);

/* val[i] = -1 means val[i] to be tuned */
int btune_tune (struct btune *btune, int *val, double *f);

int btune_cleanup (struct btune *btune);

#endif

