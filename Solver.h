//
//  Solver.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Solver_h
#define TRANSPORT_Solver_h

#include "Problem.h"
#include "Sweep.h"
#include "Mapping.h"
#include "Common.h"
#include "Region.h"

class Solver {
protected:
    vector<vector<mat> > psi_pos;
    vector<vector<mat> > psi_neg;
    vector<vec > phi;
public:
    vector<vector<mat> > get_psi_pos() {return psi_pos;};
    vector<vector<mat> > get_psi_neg() {return psi_neg;};
    vector<vec > get_phi() {return phi;};
    Solver(Input, Problem);
};


#endif
