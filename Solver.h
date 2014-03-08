//
//  Solver.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Solver_h
#define TRANSPORT_Solver_h

#include "Common.h"
#include "Problem.h"
#include "Sweep.h"
#include "Mapping.h"
#include "Fixup.h"
#include "BC.h"
#include "Adapt.h"
#include "Output.h"

void Solver (Input, Problem);

vector<vector<LDFE_reg*> > find_quad (vector<Region*>, string);

class CompareByDirs {
public:
	bool operator() (LDFE_reg* a, LDFE_reg* b) {
		vec dirs_a = a->get_dirs();
		vec dirs_b = b->get_dirs();
		if (dirs_a(0)>0) {
			return dirs_a(0) < dirs_b(0);
		}
		else {
			return dirs_a(0) > dirs_b(0);
		}
	}
};
#endif
