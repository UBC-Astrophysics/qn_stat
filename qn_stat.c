/*
 *  Copyright (C) 2015 Enterome
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "qn_stat.h"
#include <stdbool.h>
#include <string.h>

// Static functions declarations
static int cmp_double(const void*, const void*);
static double whimed(double*, uint64_t*, size_t, qn_stat_workspace*);
static double select_kth_element(double*, size_t, size_t, double*);


double qn_calc (const double* x, const size_t n) {
  double res;
  qn_stat_workspace* workspace = qn_stat_alloc(n);

  res = qn_stat(x,n,workspace);
  qn_stat_free(workspace);
  return(res);
}


qn_stat_workspace* qn_stat_alloc(const size_t n_max)
{
	qn_stat_workspace* workspace = 
		(qn_stat_workspace*)malloc(sizeof(qn_stat_workspace));

	workspace->y = (double*)malloc(n_max*sizeof(double));
	workspace->work = (double*)malloc(n_max*sizeof(double));

	workspace->left = (size_t*)malloc(n_max*sizeof(size_t));
	workspace->right = (size_t*)malloc(n_max*sizeof(size_t));
	workspace->p = (size_t*)malloc(n_max*sizeof(size_t));
	workspace->q = (size_t*)malloc(n_max*sizeof(size_t));
	workspace->weight = (uint64_t*)malloc(n_max*sizeof(uint64_t));

	workspace->a_cand = (double*)malloc(n_max*sizeof(double));
	workspace->w_cand = (double*)malloc(n_max*sizeof(double));

	workspace->b = (double*)malloc(n_max*sizeof(double));

	return workspace;
}

void qn_stat_free(qn_stat_workspace* workspace)
{
	free(workspace->y);
	free(workspace->work);

	free(workspace->left);
	free(workspace->right);
	free(workspace->p);
	free(workspace->q);
	free(workspace->weight);

	free(workspace->a_cand);
	free(workspace->w_cand);

	free(workspace->b);

	free(workspace);
}

double qn_stat(const double* x, const size_t n, qn_stat_workspace* workspace)
{
	double ret_val;

	double* y = workspace->y;
	double* work = workspace->work;

	size_t* left = workspace->left;
	size_t* right = workspace->right;
	size_t* p = workspace->p;
	size_t* q = workspace->q;
	size_t* weight = workspace->weight;

	double trial;
	bool found;

	size_t h, j, jh;
	uint64_t sump, sumq;

	double dn;
	size_t nl, nr;
	size_t k, knew;

	/* Function Body */
	h = n / 2 + 1;
	k = h * (h - 1) / 2;
	memcpy(y, x, n*sizeof(double));
	qsort(y, n, sizeof(double), cmp_double);

	for (size_t i = 0; i < n; ++i)
	{
		left[i] = n - i + 1;
		right[i] = (i <= h) ? n : n - (i - h);
	}

    nl = n * (n + 1) / 2;
    nr = n * n;
    knew = k + nl;/* = k + (n+1 \over 2) */
	found = false;

	while(!found && nr - nl > n)
	{
		j = 0;

		for (size_t i = 1; i < n; ++i)
		{
			if (left[i] <= right[i])
			{
				weight[j] = right[i] - left[i] + 1;
				jh = left[i] + weight[j] / 2;
				work[j] = y[i] - y[n - jh];
				++j;
			}
		}
		trial = whimed(work, weight, j, workspace);

		j = 0;
		for (size_t i = n; i--;)
		{
			while (j < n && (y[i] - y[n - j - 1]) < trial)
				++j;
			p[i] = j;
		}

		j = n + 1;
		for (size_t i = 0; i < n; ++i)
		{
			while ((y[i] - y[n - j + 1]) > trial)
				--j;
			q[i] = j;
		}
		sump = 0;
		sumq = 0;
		for (size_t i = 0; i < n; ++i)
		{
			sump += p[i];
			sumq += q[i] - 1;
		}

		if (knew <= sump)
		{
			for (size_t i = 0; i < n; ++i)
			{
				right[i] = p[i];
			}
			nr = sump;
		} 
		else if (knew > sumq)
		{
			for (size_t i = 0; i < n; ++i)
			{
				left[i] = q[i];
			}
			nl = sumq;
		} 
		else /* sump < knew <= sumq */
		{
			found = true;
		}
	} /* while */

	if (! found)
	{
		j = 0;
		for (size_t i = 1; i < n; ++i)
		{
			for (size_t jj = left[i]; jj <= right[i]; ++jj)
			{
				work[j] = (y[i] - y[n - jj]);
				j++;
			}/* j will be = sum_{i=2}^n (right[i] - left[i] + 1)_{+}  */
		}
		ret_val = select_kth_element(work, j, knew-nl-1, workspace->b);
	}
	else
	{
		ret_val = trial;
	}

	if (n <= 9)
	{
		if (n == 2)
		{
			dn = .399;
		}
		else if (n == 3)
		{
			dn = .994;
		}
		else if (n == 4)
		{
			dn = .512;
		}
		else if (n == 5)
		{
			dn = .844;
		}
		else if (n == 6)
		{
			dn = .611;
		}
		else if (n == 7)
		{
			dn = .857;
		}
		else if (n == 8)
		{
			dn = .669;
		}
		else 
		{
			dn = .872;
		}
	}
	else
	{
		if (n % 2 == 1)
		{
			dn = n / (n + 1.4);
		}
		else
		{
			dn = n / (n + 3.8);
		}
	}

	ret_val = dn * 2.2219 * ret_val;
	return ret_val;
}

static int cmp_double(const void *x, const void *y)
{
	const double xx = *(double*)x; 
	const double yy = *(double*)y;

	if (xx < yy) 
		return -1;
	else if (xx > yy) 
		return  1;
	else
		return 0;
}

static double whimed(double* a, uint64_t* w, size_t n, qn_stat_workspace* workspace)
{
	double* a_cand = workspace->a_cand;
	double* w_cand = workspace->w_cand;

	size_t kcand;

	uint64_t wleft, wmid, wright, w_tot, wrest;

	double trial;

	w_tot = 0;
	for (size_t i = 0; i < n; ++i)
	{
		w_tot += w[i];
	}

	wrest = 0;
	do
	{
		trial = select_kth_element(a, n, n/2, workspace->b);

		wleft = 0; wmid = 0; wright= 0;
		for (size_t i = 0; i < n; ++i)
		{
			if (a[i] < trial)
				wleft += w[i];
			else if (a[i] > trial)
				wright += w[i];
			else
				wmid += w[i];
		}

		kcand = 0;
		if (2 * (wrest + wleft) > w_tot)
		{
			for (size_t i = 0; i < n; ++i)
			{
				if (a[i] < trial)
				{
					a_cand[kcand] = a[i];
					w_cand[kcand] = w[i]; ++kcand;
				}
			}
		}
		else if (2 * (wrest + wleft + wmid) <= w_tot)
		{
			for (size_t i = 0; i < n; ++i)
			{
				if (a[i] > trial)
				{
					a_cand[kcand] = a[i];
					w_cand[kcand] = w[i]; 
					++kcand;
				}
			}
			wrest += wleft + wmid;
		}
		else
		{
			return trial;
		}

		n = kcand;
		for (size_t i = 0; i < n; ++i)
		{
			a[i] = a_cand[i];
			w[i] = w_cand[i];
		}
	} while(true);
}

static double select_kth_element(double *a, size_t n, size_t k, double* b)
{
	int64_t j, l, lr, jnc;
	double ax, buffer;

	memcpy(b, a, n*sizeof(double));

	l = 0;
	lr = n-1;

	while (l < lr) 
	{
		ax = b[k];
		jnc = l;
		j = lr;
		while (jnc <= j)
		{
			while (b[jnc] < ax)
			{
				++jnc;
			}
			while (b[j] > ax)
			{
				--j;
			}

			if (jnc <= j)
			{
				buffer = b[jnc];
				b[jnc] = b[j];
				b[j] = buffer;
				++jnc;
				--j;
			}
		}

		if (j < k)
		{
			l = jnc;
		}

		if (k < jnc)
		{
			lr = j;
		}
	}

	return b[k];
}

