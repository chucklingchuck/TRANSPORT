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
    double si_tol = input_param.get_si_tol();
    int max_si_cycles = input_param.get_si_cycles();
    double si_error;
    int si_cycles;
    Col<int> scat_order = input_param.get_scat_order();
    vec scat_xs = input_param.get_scat_xs();

    /* Initialize basic parameters */
    int num_reg = int(Regions.size());
    vector<int> num_pos(num_reg);
    vector<int> num_neg(num_reg);
    vector<mat> quad_dirs(num_reg);
    vector<mat> quad_wgts(num_reg);
    vector<int> num_cells(num_reg);
    vector<int> num_edges(num_reg);
    vector<int> num_ang(num_reg);
    string dir;
    vector<Sweep*> pos_spat_ptr(num_reg);
    vector<Sweep*> neg_spat_ptr(num_reg);
    Mapping* map_ptr;
    double leg;
    vector<vector<mat> > phi_old;
    vector<vector<mat> > phi_old_old;
    double sigma;
    double sum_top;
    double sum_bot;
    int counter;
    
    /* Calculate number of positive and negative angles in each region */
    for (int i=0; i<num_reg; i++) {
        num_pos[i] = 0;
        quad_dirs[i] = Regions[i]->get_quad_dirs();
        quad_wgts[i] = Regions[i]->get_quad_wgts();
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
        cout << quad_dirs[i] << endl;
        cout << quad_wgts[i] << endl;
    }
    
    /* Resize positive and negative psi and phi vectors to number of edges */
    psi_pos.resize(num_reg);
    psi_neg.resize(num_reg);
    phi.resize(num_reg);
    for (int i=0; i<num_reg; i++) {
        psi_pos[i].resize(num_pos[i]);
        psi_neg[i].resize(num_neg[i]);
        phi[i].resize(scat_order(i)+1);
        for (int j=0; j<num_pos[i]; j++) {
            psi_pos[i][j].resize(2, num_edges[i]);
        }
        for (int j=0; j<num_neg[i]; j++) {
            psi_neg[i][j].resize(2, num_edges[i]);
        }
        for (int j=0; j<=scat_order(i); j++){
            phi[i][j] = zeros(2, num_edges[i]-1);
        }
    }
    
    /* Assign left boundary condition */
    for (int i=0; i<num_pos[0]; i++) {
        psi_pos[0][i](0,0) = psi_left;
    }
    
    /* Assign right boundary condition */
    for (int i=0; i<num_neg[num_reg-1]; i++) {
        psi_neg[num_reg-1][i](0,num_edges[num_reg-1]-1) = psi_right;
    }
    
    /* Source iteration */
    si_error = si_tol * 2.;
    si_cycles = 0;
    do {
        
        /* Positive sweep */
        dir = "pos";
        for (int i=0; i<num_reg; i++) {
            
            /* Sweep using selected spatial discretization */
            if (spat_type == "LDFE") {
                pos_spat_ptr[i] = new LDFE_spat(Regions[i], psi_pos[i], phi[i], dir);
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
        for (int i=num_reg-1; i>=0; i--) {
            
            /* Sweep using selected spatial discretization */
            if (spat_type == "LDFE") {
                neg_spat_ptr[i] = new LDFE_spat(Regions[i], psi_neg[i], phi[i], dir);
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
        
        /* Calculate new flux moments */
        if (si_cycles>0) {
            phi_old_old = phi_old;
        }
        phi_old = phi;
        for (int i=0; i<num_reg; i++) {
            /* Rearrange quadrature directions */
            vec quad_dirs_use(num_ang[i]);
            vec quad_wgts_use(num_ang[i]);
            counter=0;
            for (int m=0; m<num_ang[i]/2; m++) {
                for (int n=0; n<2; n++) {
                    quad_dirs_use(counter)=quad_dirs[i](m,n);
                    quad_wgts_use(counter)=quad_wgts[i](m,n);
                    counter++;
                }
            }
            for (int j=0; j<=scat_order[i]; j++) {
                phi[i][j]=zeros(2,num_edges[i]-1);
                for (int k=0; k<num_edges[i]-1; k++) {
                    /* Positive angle contributions */
                    for (int r=0; r<num_pos[i]; r++) {
                        leg=legendre(double(j),quad_dirs_use(r));
                        phi[i][j](0,k)=phi[i][j](0,k)+psi_pos[i][r](1,k)*leg*quad_wgts_use(r);
                        phi[i][j](1,k)=phi[i][j](1,k)+psi_pos[i][r](0,k+1)*leg*quad_wgts_use(r);
                    }
                    /* Negative angle contributions */
                    counter=0;
                    for (int r=num_pos[i]; r<num_ang[i]; r++) {
                        leg=legendre(double(j),quad_dirs_use(r));
                        phi[i][j](0,k)=phi[i][j](0,k)+psi_neg[i][counter](0,k)*leg*quad_wgts_use(r);
                        phi[i][j](1,k)=phi[i][j](1,k)+psi_neg[i][counter](1,k+1)*leg*quad_wgts_use(r);
                        counter++;
                    }
                }
            }
        }
        
        /* Calculate new SI error */
        if (si_cycles>0) {
            sum_top=0;
            sum_bot=0;
            for (int i=0; i<num_reg; i++) {
                for (int j=0; j<num_cells[i]; j++) {
                    sum_top=sum_top+pow((phi[i][0](0,j)+phi[i][0](1,j))/2-\
                                        (phi_old[i][0](0,j)+phi_old[i][0](1,j))/2,2);
                    sum_bot=sum_bot+pow((phi_old[i][0](0,j)+phi_old[i][0](1,j))/2-\
                                        (phi_old_old[i][0](0,j)+phi_old_old[i][0](1,j))/2,2);
                }
            }
            sigma=sqrt(sum_top)/sqrt(sum_bot);
            si_error=sqrt(sum_bot)/(1.-sigma);
        }
        si_cycles++;
        //cout << si_error << endl;
        //cout << si_cycles << endl;
        
    } while (si_error > si_tol && si_cycles < max_si_cycles && scat_xs.max() > 0.);

    /* Check if SI converged */
    if (si_error > si_tol) {
        exit(EXIT_FAILURE);
    }
}
