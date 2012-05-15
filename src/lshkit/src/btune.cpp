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

#include <cstdlib>
#include <cassert>
#include <cstring>
#include "btune.h"

struct fv
{
	double	f[2];
	int	v[0];
};

#define SIZEOF_FV(N) (sizeof(struct fv) + sizeof(int) * N)
#define ALLOCA_FV(N) ((struct fv *)alloca(SIZEOF_FV(N)))

using namespace std;

int btune_init (struct btune *btune, int N, const struct btune_minmax *minmax, btune_f fun, void *param)
{
	if (db_create(&btune->db, NULL, 0) != 0) assert(0);
	if (btune->db->open(btune->db, NULL, NULL, NULL, DB_HASH, DB_CREATE, 0) != 0) assert(0);
	btune->N = N;
	btune->minmax = (struct btune_minmax *)malloc(N * sizeof(struct btune_minmax));
	memcpy(btune->minmax, minmax, N * sizeof(struct btune_minmax));
	btune->fun = fun;
	btune->param = param;
}

int btune_cleanup (struct btune *btune)
{
	btune->db->close(btune->db, 0);
	free(btune->minmax);
	return 0;
}

int btune_db_get (struct btune *btune, const int *val, struct fv *fv)
{
	DB *db = btune->db;
	DBT key, data;
	memset(&key, 0, sizeof key);
	memset(&data, 0, sizeof data);
	key.data = (void *)val;
	key.size = sizeof(int) * btune->N;
	data.data = fv;
	data.ulen = SIZEOF_FV(btune->N);
	data.flags = DB_DBT_USERMEM;
	return db->get(db, NULL, &key, &data, 0);
}

int btune_db_put (struct btune *btune, const int *val, const struct fv *fv)
{
	DB *db = btune->db;
	DBT key, data;
	memset(&key, 0, sizeof key);
	memset(&data, 0, sizeof data);
	key.data = (void *)val;
	key.size = sizeof(int) * btune->N;
	data.data = (void *)fv;
	data.size = SIZEOF_FV(btune->N);
	if (db->put(db, NULL, &key, &data, 0) != 0) assert(0);
	return 0;
}

/* 1 if succeed, 0 if fail */
int btune_help (struct btune *btune, int *val, struct fv *fv, int n, int d)
{
	if (btune_db_get(btune, val, fv) == 0)
	{
		return fv->f[1] >= 0;
	}
	else if (n > d)
	{
		btune->fun(fv->f, val, btune->param);
		memcpy(fv->v, val, sizeof(int) * btune->N);
		btune_db_put(btune, val, fv);
		return fv->f[1] >= 0.0;
	}
	else if (val[n] != -1) return btune_help(btune, val, fv, n+1, d);
	else	/* tune val[n] */
	{
		int min = btune->minmax[n].min;
		int max = btune->minmax[n].max;
		int t;
		double vt;

		/* make sure max works */
		val[n] = max;
		if (btune_help(btune, val, fv, n+1, d) == 0)
		{
			val[n] = -1;
			btune_db_put(btune, val, fv);
			return 0;
		}

		/* locate the mininal min that works */
		t = min;
		while (t <= max)
		{
			val[n] = t;
			if (btune_help(btune, val, fv, n+1, d)) break;
			t *= 2;
			if (t == 0) t = 1;
		}
		if (t > min)
		{
			int r, s;
			/* s doen't work, t works, s < t */
			s = t / 2;
			if (t > max) t = max;
			while (s + 1 < t)
			{
				r = (s + t) / 2;
				val[n] = r;
				if (btune_help(btune, val, fv, n+1, d)) t = r;
				else s = r;
			}
			min = t;
		}

		if (min >= max)
		{
			val[n] = min;
			btune_help(btune, val, fv, n+1, d);
			val[n] = -1;
			btune_db_put(btune, val, fv);
			return 1;
		}

		/* both min, max works, locate optimal */
		val[n] = t = min + 1;
		if (btune_help(btune, val, fv, n+1, d) == 0) assert(0);
		vt = fv->f[0];

		val[n] = min;
		if (btune_help(btune, val, fv, n+1, d) == 0) assert(0);
		if (fv->f[0] < vt)
		{
			val[n] = -1;
			btune_db_put(btune, val, fv);
			return 1;
		}

		min = t;
		for (;;)
		{
			int s;
			double vs;
			if (min >= max) break;
			t = (min + max + 1) / 2;
			val[n] = t;
			if (btune_help(btune, val, fv, n+1, d) == 0) assert(0);
			vt = fv->f[0];

			s = t - 1;
			val[n] = s;
			if (btune_help(btune, val, fv, n+1, d) == 0) assert(0);
			vs = fv->f[0];
			if (vs < vt)
			{
				max = s;
			}
			else
			{
				min = t;
			}
		}
		assert(min == max);
		val[n] = min;
		if (btune_help(btune, val, fv, n+1, d) == 0) assert(0);
		val[n] = -1;
		btune_db_put(btune, val, fv);
		return 1;
	}
}

int btune_tune (struct btune *btune, int *val, double *f)
{
	struct fv *fv;
	int i = btune->N - 1;
	int r;
	while ((i >= 0) && (val[i] != -1)) i--;
	if (i < 0) return -1;
	fv = ALLOCA_FV(btune->N);
	r = btune_help(btune, val, fv, 0, i);
	if (r == 1)
	{
		f[0] = fv->f[0];
		f[1] = fv->f[1];
		memcpy(val, fv->v, sizeof(int) * btune->N);
	}
	return r;
}

