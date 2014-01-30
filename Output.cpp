//
//  Output.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/25/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Output.h"

void Output(Problem problem_setup, Solver solution) {
    
    /* Define basic parameters used */
    vector<Region*> Regions = problem_setup.get_Regions();
    int num_reg = int(Regions.size());
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
    
    /* Write quadratures of each region to separate HDF5 files */
    for (int i=0; i<num_reg; i++) {
        mat quad_dirs = Regions[i]->get_quad_dirs();
        mat quad_wgts = Regions[i]->get_quad_wgts();
        snprintf(buffer, sizeof(char) * 32, "dirs_reg%i", i);
        quad_dirs.save(buffer, hdf5_binary);
        quad_dirs.clear();
        snprintf(buffer, sizeof(char) * 32, "wgts_reg%i", i);
        quad_wgts.save(buffer, hdf5_binary);
        quad_wgts.clear();
    }
    
    /* Write other problem parameters to separate HDF5 files */
    vec reg_size(num_reg);
    Col<int> num_cells(num_reg);
    vec abs_xs(num_reg);
    vec ext_source(num_reg);
    for (int i=0; i<num_reg; i++) {
        reg_size(i) = Regions[i]->get_reg_size();
        num_cells(i) = Regions[i]->get_num_cells();
        abs_xs(i) = Regions[i]->get_abs_xs();
        ext_source(i) = Regions[i]->get_ext_source();
    }
    reg_size.save("reg_size", hdf5_binary);
    num_cells.save("num_cells", hdf5_binary);
    abs_xs.save("abs_xs",hdf5_binary);
    ext_source.save("ext_source",hdf5_binary);
    
}