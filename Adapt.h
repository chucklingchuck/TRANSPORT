/*
 * Adapt.h
 *
 *  Created on: Feb 28, 2014
 *      Author: cheuklau
 */

#ifndef ADAPT_H_
#define ADAPT_H_

#include "Common.h"
#include "Region.h"
#include "Solver.h"

void adapt(vector<vector<LDFE_reg*> >, vector<vector<LDFE_reg*> >, Col<int>, string, double, double, vector<Region*>, bool&, double, double);

void apply_adapt_ind(vector<LDFE_reg*>, vector<vec>&,  int, string, double);

void coarsen(vector<Region*>, string);

#endif /* ADAPT_H_ */
