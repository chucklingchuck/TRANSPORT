//
//  Solver.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Solver.h"

Solver::Solver (Input input_param, Problem problem_setup) {
    
    /* Get region information and boundary conditions */
    vector<Region*> Regions = problem_setup.get_Regions();
    double psi_left = input_param.get_psi_left();
    double psi_right = input_param.get_psi_right();
    string spat_type = input_param.get_spat_type();
    
    /* Initialize basic parameters */
    int num_reg = int(Regions.size());
    vector<int> num_pos(num_reg);
    vector<int> num_neg(num_reg);
    vector<mat> quad_dirs(num_reg);
    vector<int> num_cells(num_reg);
    vector<int> num_edges(num_reg);
    vector<int> num_ang(num_reg);
    
    /* Calculate number of positive and negative angles in each region */
    for (int i=0; i<num_reg; i++) {
        num_pos[i] = 0;
        quad_dirs[i] = Regions[i]->get_quad_dirs();
        num_cells[i] = Regions[i]->get_num_cells();
        num_edges[i] = num_cells[i] + 1;
        num_ang[i] = quad_dirs[i].n_rows * 2;
        for (int j=0; j<num_ang[i]/2; j++) {
            if (quad_dirs[i](j,1)>0) {
                num_pos[i]++;
            }
        }
        num_pos[i] = num_pos[i]*2;
        num_neg[i] = num_ang[i]-num_pos[i];
    }
    
    /* Resize positive and negative psi and phi vectors to number of edges */
    psi_pos.resize(num_reg);
    psi_neg.resize(num_reg);
    phi.resize(num_reg);
    for (int i=0; i<num_reg; i++) {
        psi_pos[i].resize(num_pos[i]);
        psi_neg[i].resize(num_neg[i]);
        for (int j=0; j<num_pos[i]; j++) {
            psi_pos[i][j].resize(2, num_edges[i]);
        }
        for (int j=0; j<num_neg[i]; j++) {
            psi_neg[i][j].resize(2, num_edges[i]);
        }
        phi[i].resize(num_edges[i]);
    }
    
    /* Assign left boundary condition */
    for (int i=0; i<num_pos[0]; i++) {
        psi_pos[0][i](0,0) = psi_left;
    }
    
    /* Assign right boundary condition */
    for (int i=0; i<num_neg[num_reg-1]; i++) {
        psi_neg[num_reg-1][i](0,num_edges[num_reg-1]-1) = psi_right;
    }
    
    /* Positive sweep */
    string dir = "pos";
    vector<Sweep*> pos_spat_ptr(num_reg);
    Mapping* map_ptr;
    for (int i=0; i<num_reg; i++) {
        
        /* Sweep using selected spatial discretization */
        if (spat_type == "LDFE") {
            pos_spat_ptr[i] = new LDFE_spat(Regions[i], psi_pos[i], dir);
        }
        else {
            cout << "Invalid spatial discretization type. Exiting program." << endl;
            exit (EXIT_FAILURE);
        }
        
        /* Store positive psi solution */
        psi_pos[i] = pos_spat_ptr[i]->get_region_psi();
        
        /* Map solution to downstream region */
        if (i<num_reg-1) {
            
            /* Calculate current and downstream quadrature order */
            int current_order = Regions[i]->get_order();
            int downstream_order = Regions[i+1]->get_order();
            
            /* One-to-one mapping */
            if (current_order == downstream_order) {
                map_ptr = new one_to_one(psi_pos[i], dir);
            }
            
            /* Fine to coarse mapping */
            else if (current_order > downstream_order) {
                
            }
            
            /* Coarse to fine mapping */
            else if (current_order < downstream_order) {
            
            }
            
            /* Mapping error */
            else {
                cout << "Mapping error. Exiting program." << endl;
                exit(EXIT_FAILURE);
            }
            
            /* Store mapped solution in next region */
            vec mapped_sol = map_ptr->get_psi_to();
            for (int j=0; j<num_pos[i+1]; j++) {
                psi_pos[i+1][j](0,0) = mapped_sol(j);
            }
        }
    }
    
    /* Negative sweep */
    dir = "neg";
    vector<Sweep*> neg_spat_ptr(num_reg);
    for (int i=num_reg-1; i>=0; i--) {
        
        /* Sweep using selected spatial discretization */
        if (spat_type == "LDFE") {
            neg_spat_ptr[i] = new LDFE_spat(Regions[i], psi_neg[i], dir);
        }
        else {
            cout << "Invalid spatial discretization type. Exiting program." << endl;
            exit (EXIT_FAILURE);
        }
        
        /* Store negative psi solution */
        psi_neg[i] = neg_spat_ptr[i]->get_region_psi();
        
        /* Map solution to downstream region */
        if (i>0) {
            
            /* Calculate current and downstream quadrature order */
            int current_order = Regions[i]->get_order();
            int downstream_order = Regions[i-1]->get_order();
            
            /* One-to-one mapping */
            if (current_order == downstream_order) {
                map_ptr = new one_to_one(psi_neg[i], dir);
            }
            
            /* Fine to coarse mapping */
            else if (current_order > downstream_order) {
                
            }
            
            /* Coarse to fine mapping */
            else if (current_order < downstream_order) {
                
            }
            
            /* Mapping error */
            else {
                exit(EXIT_FAILURE);
            }
            
            /* Store mapped solution in next region */
            vec mapped_sol = map_ptr->get_psi_to();
            for (int j=0; j<num_neg[i-1]; j++) {
                psi_neg[i-1][j](0,num_edges[i-1]-1) = mapped_sol(j);
            }
        }
    }
}
