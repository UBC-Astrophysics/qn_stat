#ifndef QN_STAT_H
#define QN_STAT_H

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

#include <stdlib.h>
#include <stdint.h>

typedef struct qn_stat_workspace
{
	double* y;
	double* work;
	size_t* left;
	size_t* right;
	size_t* p;
	size_t* q;
	uint64_t* weight;

	double* a_cand;
	double* w_cand;

	double* b;
} qn_stat_workspace;

 
qn_stat_workspace* qn_stat_alloc(const size_t n_max);
void qn_stat_free(qn_stat_workspace* workspace);
double qn_stat(const double* x, const size_t n, qn_stat_workspace* workspace);

#endif // QN_STAT_H
