//
//  Solver.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Solver.h"

void Solver (Input input_param, Problem problem_setup) {
    
    /* Retrieve needed input parameters */
    int max_si_cycles       = input_param.get_si_cycles();
    double psi_left         = input_param.get_psi_left();
    double psi_right        = input_param.get_psi_right();
    double si_tol           = input_param.get_si_tol();
    Col<int> num_cells      = input_param.get_num_cells();
    Col<int> scat_order     = input_param.get_scat_order();
    vector<Region*> Regions = problem_setup.get_Regions();
    vec scat_xs             = input_param.get_scat_xs();
    vec ext_source          = input_param.get_ext_source();
    string map_type         = input_param.get_map_type();
    string fixup_type       = input_param.get_fixup_type();
    string spat_type        = input_param.get_spat_type();
    
    /* Calculate basic parameters */
    int num_reg = int(Regions.size());
    
    /* Initialize parameters */
    int si_cycles;
    double si_error;
    double sigma;
    double sum_top;
    double sum_bot;
    double leg;
    double abs_xs_reg;
    bool reached_bottom;
    string dir;
    LDFE_quad* quad_region;
    LDFE_reg* quad_selected;
    LDFE_quad* quad_ptr;
    vec dirs_use;
    vec wgts_use;
    vec tot_source_lft(num_reg - 1);
    vec tot_source_rgt(num_reg - 1);
    vec scat_xs_reg;
    vector<vector<mat> > phi_old;
    vector<vector<mat> > phi_old_old;
    vector<vector<mat> > phi;
    vector<LDFE_reg*> quad_ind;
    vector<LDFE_reg*> quad_children;
    vector<LDFE_reg*> quad_children_tmp_ind;
    vector<LDFE_reg*> quad_children_tmp_all;
    
    /* Resize phi vector */
    phi.resize(num_reg);
    for (int i=0; i<num_reg; i++) {
        phi[i].resize(scat_order(i)+1);
        for (int j=0; j<=scat_order(i); j++){
            phi[i][j] = zeros(2, num_cells(i));
        }
    }
    
    /* Source iteration */
    si_error = si_tol * 2.;
    si_cycles = 0;
    do {
        
        /* Calculate the left and right total source at spatial interfaces if using FIXUP1 */
        if (fixup_type == "FIXUP1") {
            // Go through each region interface
            for (int i=0; i<num_reg-1; i++) {
                tot_source_lft(i) = 0;                   // Initialize left scattering source for current interface
                scat_xs_reg = Regions[i]->get_scat_xs(); // Scattering cross section for left region
                abs_xs_reg = Regions[i]->get_abs_xs();   // Absorption cross section for left region
                // Add the contributions from each moment
                for (int s=0; s<=scat_order[i]; s++) {
                    leg = legendre(double(s), 0.0);
                    tot_source_lft(i)=tot_source_lft(i)+((2.*s+1.)/2.)*scat_xs_reg(s)*phi[i][s](0,num_cells(i)-1)*leg;
                }
                // Add the external source
                tot_source_lft(i) = (tot_source_lft(i)+ext_source(i))/(abs_xs_reg+scat_xs_reg(0));
                // Calculate right scattering source
                tot_source_rgt(i)=0;
                scat_xs_reg = Regions[i+1]->get_scat_xs();
                abs_xs_reg = Regions[i+1]->get_abs_xs();
                // Add the contributions from each moment
                for (int s=0; s<=scat_order[i+1]; s++) {
                    leg = legendre(double(s), 0.0);
                    tot_source_rgt(i)=tot_source_rgt(i)+((2.*s+1.)/2.)*scat_xs_reg(s)*phi[i+1][s](0,0)*leg;
                }
                // Add the external source
                tot_source_rgt(i) = (tot_source_rgt(i)+ext_source(i+1))/(abs_xs_reg+scat_xs_reg(0));
            }
        }
        /* Apply left boundary condition */
        quad_region   = Regions[0]->get_quad();        // Pointer to left spatial region quadrature
        quad_ind      = quad_region->get_quad();   // Vector of positive and negative LDFE regions for left spatial region
        quad_selected = quad_ind[0];                   // Positive LDFE regions root node for left region
        if (quad_selected->get_has_children() == true) { // Check if the coarsest quadrature was chosen
            quad_children = quad_selected->get_children(); // Root node's first level children
            // Traverse the tree
            reached_bottom = false;
            while (reached_bottom == false) {
                // Go through the children nodes in the current level
                for (int i=0; i<quad_children.size(); i++) {
                    // If the child has children add them to a temporary children vector
                    if (quad_children[i]->get_has_children() == true) {
                        quad_children_tmp_ind = quad_children[i]->get_children();
                        quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                        quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                    }
                    // If the child does not have children resize its psi and apply the left boundary condition
                    else {
                        quad_children[i]->psi.resize(2);
                        for (int j=0; j<2; j++) {
                            quad_children[i]->psi[j].resize(2, num_cells[0]+1);
                            quad_children[i]->psi[j](0,0)=psi_left;
                        }
                    }
                }
                // If there are no children in the temporary children vector the bottom is reached
                if (quad_children_tmp_all.size() == 0) {
                    reached_bottom = true;
                }
                // Set the children vector as the temporary children vector for the next level calculation
                quad_children = quad_children_tmp_all;
                quad_children_tmp_all.clear();
                quad_children_tmp_ind.clear();
            }
        }
        // Directly apply the left boundary condition for the coarsest quadrature
        else {
            quad_selected->psi.resize(2);
            for (int j=0; j<2; j++) {
                quad_selected->psi[j].resize(2, num_cells[0]+1);
                quad_selected->psi[j](0,0)=psi_left;
            }
        }
        
        /* Positive sweep */
        dir = "pos";
        for (int i=0; i<num_reg; i++) {
            // Sweep using selected spatial discretization
            if (spat_type == "LDFE") {
                LDFE_spat(Regions[i], phi[i], dir);
            }
            else {
                cout << "Invalid spatial discretization type. Exiting program." << endl;
                exit (EXIT_FAILURE);
            }
            // Map solution to downstream region
            if (i<num_reg-1) {
                // mapping1 mapping
                if (map_type == "MAPPING1") {
                    mapping1(Regions[i], Regions[i+1], dir);
                }
                else if (map_type == "MAPPING2") {
                    mapping2(Regions[i], Regions[i+1], dir);
                }
                else {
                    cout << "Invalid mapping type. Exiting program." << endl;
                    exit (EXIT_FAILURE);
                }
                // Apply fixup
                if (fixup_type == "FIXUP1") {
                    fixup1(Regions[i], Regions[i+1], dir, tot_source_rgt(i));
                }
            }
            
        }
        
        /* Apply right boundary condition */
        quad_region   = Regions[num_reg-1]->get_quad(); // Pointer to right spatial region quadrature
        quad_ind      = quad_region->get_quad();    // Vector of positive and negative LDFE regions for right spatial region
        quad_selected = quad_ind[1];                    // Negative LDFE regions root node for right region
        if (quad_selected->get_has_children() == true) {
            quad_children = quad_selected->get_children(); // Root node's first level children
            // Traverse the tree
            reached_bottom = false;
            while (reached_bottom == false) {
                // Go through the children nodes in the current level
                for (int i=0; i<quad_children.size(); i++) {
                    // If the child has children add them to a temporary children vector
                    if (quad_children[i]->get_has_children() == true) {
                        quad_children_tmp_ind = quad_children[i]->get_children();
                        quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                        quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                    }
                    // If the child does not have children resize its psi and apply the right boundary condition
                    else {
                        quad_children[i]->psi.resize(2);
                        for (int j=0; j<2; j++) {
                            quad_children[i]->psi[j].resize(2, num_cells[num_reg-1]+1);
                            quad_children[i]->psi[j](0,num_cells[num_reg-1])=psi_right;
                        }
                    }
                }
                // If there are no children in the temporary children vector the bottom is reached
                if (quad_children_tmp_all.size() == 0) {
                    reached_bottom = true;
                }
                // Set the children vector as the temporary children vector for the next level calculation
                quad_children = quad_children_tmp_all;
                quad_children_tmp_all.clear();
                quad_children_tmp_ind.clear();
            }
        }
        // Directly apply the right boundary condition for the coarsest quadrature
        else {
            quad_selected->psi.resize(2);
            for (int j=0; j<2; j++) {
                quad_selected->psi[j].resize(2, num_cells[num_reg-1]+1);
                quad_selected->psi[j](0,num_cells[num_reg-1])=psi_right;
            }
        }
        
        /* Negative sweep */
        dir = "neg";
        for (int i=num_reg-1; i>=0; i--) {
            // Sweep using selected spatial discretization
            if (spat_type == "LDFE") {
                LDFE_spat(Regions[i], phi[i], dir);
            }
            else {
                cout << "Invalid spatial discretization type. Exiting program." << endl;
                exit (EXIT_FAILURE);
            }
            // Map solution to downstream region
            if (i>0) {
                if (map_type == "MAPPING1") {
                    mapping1(Regions[i], Regions[i-1], dir);
                }
                else if (map_type == "MAPPING2") {
                    // Call mapping
                    mapping2(Regions[i], Regions[i-1], dir);
                }
                else {
                    cout << "Invalid mapping type. Exiting program." << endl;
                    exit (EXIT_FAILURE);
                }
                // Apply fixup
                if (fixup_type == "FIXUP1") {
                    fixup1(Regions[i], Regions[i-1], dir, tot_source_rgt(i-1));
                }
            }
        }
        
        /* Update old old flux moments */
        if (si_cycles>0) {
            phi_old_old = phi_old;
        }
        
        /* Calculate new flux moments and SI error */
        if (scat_xs.max()>0) {
            phi_old = phi;
            // Go through each region
            for (int i=0; i<num_reg; i++) {
                // Go through the flux moments for each region
                for (int j=0; j<=scat_order[i]; j++) {
                    phi[i][j]=zeros(2,num_cells(i));     // Resize each flux moment to number of cells
                    // Go through each cell
                    for (int k=0; k<num_cells(i); k++) {
                        // Get to the root node of interest
                        quad_region = Regions[i]->get_quad();   // This is the current region quadrature pointer
                        quad_ind = quad_region->get_quad(); // This is the vector of pos and neg LDFE regions
                        // Go through pos and neg LDFE regions
                        for (int r=0; r<2; r++) {
                            quad_selected = quad_ind[r];                   // This is the pos or neg root LDFE node
                            if (quad_selected->get_has_children() == true) {
                                quad_children = quad_selected->get_children(); // This is the root node's 1st level children
                                reached_bottom = false;
                                while (reached_bottom == false) {
                                    for (int t=0; t<quad_children.size(); t++) { // Go through all the children nodes of current level
                                        // If the child has children add them to a temporary children vector
                                        if (quad_children[t]->get_has_children() == true) {
                                            quad_children_tmp_ind = quad_children[t]->get_children();
                                            quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                                            quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                                        }
                                        // If child does not have children add their contribution to the flux moments
                                        else {
                                            dirs_use = quad_children[t]->get_dirs();
                                            wgts_use = quad_children[t]->get_wgts();
                                            for (int s=0; s<2; s++) {
                                                leg=legendre(double(j),dirs_use(s));
                                                if (r==0){
                                                    phi[i][j](0,k)=phi[i][j](0,k)+quad_children[t]->psi[s](1,k)*leg*wgts_use(s);
                                                    phi[i][j](1,k)=phi[i][j](1,k)+quad_children[t]->psi[s](0,k+1)*leg*wgts_use(s);
                                                }
                                                else if (r==1) {
                                                    phi[i][j](0,k)=phi[i][j](0,k)+quad_children[t]->psi[s](0,k)*leg*wgts_use(s);
                                                    phi[i][j](1,k)=phi[i][j](1,k)+quad_children[t]->psi[s](1,k+1)*leg*wgts_use(s);
                                                    
                                                }
                                            }
                                        }
                                    }
                                    // If there are no children in the temporary children vector the bottom is reached
                                    if (quad_children_tmp_all.size() == 0) {
                                        reached_bottom = true;
                                    }
                                    // Set the children vector as the temporary children vector for the next level calculation
                                    quad_children = quad_children_tmp_all;
                                    quad_children_tmp_all.clear();
                                    quad_children_tmp_ind.clear();
                                }
                            }
                            // Add the contribution to the flux moment directly for the coarsest quadrature
                            else {
                                dirs_use = quad_selected->get_dirs();
                                wgts_use = quad_selected->get_wgts();
                                for (int s=0; s<2; s++) {
                                    leg=legendre(double(j),dirs_use(s));
                                    if (r==0){
                                        phi[i][j](0,k)=phi[i][j](0,k)+quad_selected->psi[s](1,k)*leg*wgts_use(s);
                                        phi[i][j](1,k)=phi[i][j](1,k)+quad_selected->psi[s](0,k+1)*leg*wgts_use(s);
                                    }
                                    else if (r==1) {
                                        phi[i][j](0,k)=phi[i][j](0,k)+quad_selected->psi[s](0,k)*leg*wgts_use(s);
                                        phi[i][j](1,k)=phi[i][j](1,k)+quad_selected->psi[s](1,k+1)*leg*wgts_use(s);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Calculate new SI error
            if (si_cycles>0) {
                sum_top=0;
                sum_bot=0;
                for (int i=0; i<num_reg; i++) {
                    for (int j=0; j<num_cells(i); j++) {
                        sum_top=sum_top+pow((phi[i][0](0,j)+phi[i][0](1,j))/2-\
                                            (phi_old[i][0](0,j)+phi_old[i][0](1,j))/2,2);
                        sum_bot=sum_bot+pow((phi_old[i][0](0,j)+phi_old[i][0](1,j))/2-\
                                            (phi_old_old[i][0](0,j)+phi_old_old[i][0](1,j))/2,2);
                    }
                }
                sigma=sqrt(sum_top)/sqrt(sum_bot);
                si_error=sqrt(sum_bot)/(1.-sigma);
            }
            cout << "current SI cycle: " << si_cycles << " current SI error: "<< si_error << endl;
            si_cycles++;
        }
        
    } while (abs(si_error) > si_tol && si_cycles < max_si_cycles && scat_xs.max() > 0.);
    
    /* Check if SI successfully converged */
    if (si_error > si_tol  && scat_xs.max() > 0.) {
        cout << "SI iterations did not converge in the requested iteration limit." << endl;
        exit(EXIT_FAILURE);
    }
}
