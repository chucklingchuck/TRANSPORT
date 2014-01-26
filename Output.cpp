//
//  Output.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/25/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Output.h"

void Output(Input input_param, Solver solution) {
    
    /* Define basic parameters used */
    int num_reg = input_param.get_num_reg();
    vector<vector<mat> > psi_pos = solution.get_psi_pos();
    vector<vector<mat> > psi_neg = solution.get_psi_neg();
    int num_pos;
    int num_neg;
    int num_edges;
    int num_total;
    char buffer[32];
    
    /* Write the angular flux of each region to separate HDF5 files */
    for (int i=0; i<num_reg; i++) {
        num_pos = int(psi_pos[i].size());
        num_neg = int(psi_neg[i].size());
        num_total = num_pos+num_neg;
        num_edges = int(psi_pos[i][0].n_cols);
        cube psi_store(2, num_edges, num_total);
        for (int j=0; j<num_pos; j++) {
            psi_store.slice(j) = psi_pos[i][j];
        }
        for (int j=num_pos; j<num_total; j++) {
            psi_store.slice(j) = psi_neg[i][j-num_pos];
        }
        snprintf(buffer, sizeof(char) * 32, "psi_reg%i", i);
        psi_store.save(buffer, hdf5_binary);
        psi_store.clear();
    }
    
}