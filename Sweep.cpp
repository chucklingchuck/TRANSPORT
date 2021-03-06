//
//  Sweep.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Sweep.h"

LDFE_spat::LDFE_spat(Region* current_reg, vector<mat> current_psi, string current_dir) {
    
    /* Define and initialize basic parameters */
    double abs_xs  = current_reg->get_abs_xs();
    double scat_xs = current_reg->get_scat_xs();
    double tot_xs = abs_xs + scat_xs;
    double reg_size = current_reg->get_reg_size();
    int num_cells = current_reg->get_num_cells();
    int num_edges = num_cells+1;
    double mesh_size = reg_size/num_cells;
    int num_angs = int(current_psi.size());
    mat quad_dirs = current_reg->get_quad_dirs();
    double tau;
    
    /* Positive sweep */
    if (current_dir == "pos") {
        
        /* Rearrange quadrature directions */
        vec quad_dirs_use(num_angs);
        int counter=0;
        for (int i=0; i<num_angs/2; i++) {
            for (int j=0; j<2; j++) {
                quad_dirs_use(counter)=quad_dirs(i,j);
                counter++;
            }
        }
        
        /* Sweep */
        for (int i=0; i<num_angs; i++) {
            tau = tot_xs*mesh_size/quad_dirs_use(i);
            for (int j=1; j<num_edges; j++) {
                current_psi[i](0,j)=current_psi[i](0,j-1)*(6.-2.*tau)/(6.+4.*tau+pow(tau,2.));
                current_psi[i](1,j-1)=current_psi[i](0,j-1)*(6.+4.*tau)/(6.+4.*tau+pow(tau,2.));
            }
        }
        
        /* Store calculated psi */
        region_psi = current_psi;
    }
    
    /* Negative sweep */
    else if (current_dir == "neg") {
        
        /* Calculate number of positive LDFE ranges */
        int num_pos_reg = quad_dirs.n_rows - (num_angs/2);
        
        /* Rearrange quadrature directions */
        vec quad_dirs_use(num_angs);
        int counter=0;
        for (int i=0; i<num_angs/2; i++) {
            for (int j=0; j<2; j++) {
                quad_dirs_use(counter)=quad_dirs(num_pos_reg+i,j);
                counter++;
            }
        }
        
        /* Sweep */
        for (int i=0; i<num_angs; i++) {
            tau = tot_xs * mesh_size / quad_dirs_use(i);
            for (int j=num_edges-2; j>=0; j--) {
                current_psi[i](0,j)=current_psi[i](0,j+1)*(6.+2.*tau)/(6.-4.*tau+pow(tau,2.));
                current_psi[i](1,j+1)=current_psi[i](0,j+1)*(6.-4.*tau)/(6.-4.*tau+pow(tau,2.));
            }
        }
        
        /* Store calculated psi */
        region_psi = current_psi;
    }
}