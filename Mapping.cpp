//
//  Mapping.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Mapping.h"

one_to_one::one_to_one (vector<mat> psi_from, string dir) {

    /* Resize mapped solution vector */
    int num_dirs = int(psi_from.size());
    int num_edges = psi_from[0].n_cols;
    psi_to.resize(num_dirs);
    
    /* Positive direction mapping */
    if (dir == "pos") {
        for (int i=0; i<num_dirs; i++) {
            psi_to(i) = psi_from[i](0, num_edges-1);
        }
    }
    
    /* Negative direction mapping */
    else if(dir == "neg") {
        for (int i=0; i<num_dirs; i++) {
            psi_to(i) = psi_from[i](0, 0);
        }
    }
}