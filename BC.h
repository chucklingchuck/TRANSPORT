/*
 * BC.h
 *
 *  Created on: Feb 27, 2014
 *      Author: cheuklau
 */

#ifndef BC_H_
#define BC_H_

#include "Common.h"
#include "Region.h"

void apply_bc(string, vector<vector<LDFE_reg*> >, vector<vector<LDFE_reg*> >, double, double, Col<int>, int, vector<Region*>, vec, vector<vec>&);

#endif /* BC_H_ */
