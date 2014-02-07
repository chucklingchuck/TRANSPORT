//
//  Sweep.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Sweep.h"

void LDFE_spat(Region* current_reg, vector<mat> phi, string current_dir) {
    
    /* Retrieve needed input parameters */
    int num_cells          = current_reg->get_num_cells();
    int scat_order         = current_reg->get_scat_order();
    double abs_xs          = current_reg->get_abs_xs();
    double reg_size        = current_reg->get_reg_size();
    double ext_source      = current_reg->get_ext_source();
    vec scat_xs            = current_reg->get_scat_xs();
    LDFE_quad* quad_region = current_reg->get_quad();

    /* Calculate basic parameters */
    double tot_xs    = abs_xs+scat_xs(0);
    double mesh_size = reg_size/num_cells;
    int num_edges    = num_cells+1;
    
    /* Initialize needed variables */
    double tau;
    double leg;
    double scat_source_lft;
    double scat_source_rgt;
    bool reached_bottom;
    vec dirs_use;
    vec wgts_use;
    vector<LDFE_reg*> quad_ind;
    vector<LDFE_reg*> quad_children;
    vector<LDFE_reg*> quad_children_tmp_ind;
    vector<LDFE_reg*> quad_children_tmp_all;
    LDFE_reg* quad_selected;

    /* Positive sweep */
    if (current_dir == "pos") {
        quad_ind      = quad_region->get_quad();   // Vector of positive and negative LDFE regions for current region
        quad_selected = quad_ind[0];                   // Positive LDFE regions root node for current region
        if (quad_selected->get_has_children()==true) {
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
                    // If the child does not have children then sweep
                    else {
                        dirs_use = quad_children[i]->get_dirs();
                        wgts_use = quad_children[i]->get_wgts();
                        // Go through each direction for each LDFE region
                        for (int j=0; j<2; j++) {
                            tau = tot_xs*mesh_size/dirs_use(j);
                            // Go through each spatial edge
                            for (int k=1; k<num_edges; k++) {
                                // Calculate scattering source
                                scat_source_lft=0;
                                scat_source_rgt=0;
                                for (int r=0; r<=scat_order; r++) {
                                    leg = legendre(double(r), dirs_use(j));
                                    scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k-1)*leg;
                                    scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k-1)*leg;
                                }
                                // Downstream angular flux
                                quad_children[i]->psi[j](0,k)=\
                                quad_children[i]->psi[j](0,k-1)*(6.-2.*tau)/(6.+4.*tau+pow(tau,2.))+\
                                ((scat_source_lft+ext_source)/tot_xs)*(3.*tau/(6.+4.*tau+pow(tau,2.)))+\
                                ((scat_source_rgt+ext_source)/tot_xs)*(3.*tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.));
                                // Upstream angular flux
                                quad_children[i]->psi[j](1,k-1)=\
                                quad_children[i]->psi[j](0,k-1)*(6.+4.*tau)/(6.+4.*tau+pow(tau,2.))+\
                                ((scat_source_lft+ext_source)/tot_xs)*(tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.))-\
                                ((scat_source_rgt+ext_source)/tot_xs)*(tau/(6.+4.*tau+pow(tau,2.)));
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
        // Directly sweep the directions for the coarsest quadrature
        else {
            dirs_use = quad_selected->get_dirs();
            wgts_use = quad_selected->get_wgts();
            // Go through each direction for the LDFE region
            for (int j=0; j<2; j++) {
                tau = tot_xs*mesh_size/dirs_use(j);
                // Go through each spatial edge
                for (int k=1; k<num_edges; k++) {
                    // Calculate scattering source
                    scat_source_lft=0;
                    scat_source_rgt=0;
                    for (int r=0; r<=scat_order; r++) {
                        leg = legendre(double(r), dirs_use(j));
                        scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k-1)*leg;
                        scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k-1)*leg;
                    }
                    // Downstream angular flux
                    quad_selected->psi[j](0,k)=\
                    quad_selected->psi[j](0,k-1)*(6.-2.*tau)/(6.+4.*tau+pow(tau,2.))+\
                    ((scat_source_lft+ext_source)/tot_xs)*(3.*tau/(6.+4.*tau+pow(tau,2.)))+\
                    ((scat_source_rgt+ext_source)/tot_xs)*(3.*tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.));
                    // Upstream angular flux
                    quad_selected->psi[j](1,k-1)=\
                    quad_selected->psi[j](0,k-1)*(6.+4.*tau)/(6.+4.*tau+pow(tau,2.))+\
                    ((scat_source_lft+ext_source)/tot_xs)*(tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.))-\
                    ((scat_source_rgt+ext_source)/tot_xs)*(tau/(6.+4.*tau+pow(tau,2.)));
                }
            }
        }
    }
    
    /* Negative sweep */
    else if (current_dir == "neg") {
        quad_ind      = quad_region->get_quad();   // Vector of positive and negative LDFE regions for current region
        quad_selected = quad_ind[1];                   // Negative LDFE regions root node for current region
        if (quad_selected->get_has_children()==true) {
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
                    // If the child does not have children then sweep
                    else {
                        dirs_use = quad_children[i]->get_dirs();
                        wgts_use = quad_children[i]->get_wgts();
                        // Go through each direction for each LDFE region
                        for (int j=0; j<2; j++) {
                            tau = tot_xs*mesh_size/dirs_use(j);
                            // Go through each spatial edge
                            for (int k=num_edges-2; k>=0; k--) {
                                // Calculate scattering source
                                scat_source_lft=0;
                                scat_source_rgt=0;
                                for (int r=0; r<=scat_order; r++) {
                                    leg = legendre(double(r), dirs_use(j));
                                    scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k)*leg;
                                    scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k)*leg;
                                }
                                // Downstream angular flux
                                quad_children[i]->psi[j](0,k)=\
                                quad_children[i]->psi[j](0,k+1)*(6.+2.*tau)/(6.-4.*tau+pow(tau,2.))+\
                                ((scat_source_lft+ext_source)/tot_xs)*(-3.*tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.))-\
                                ((scat_source_rgt+ext_source)/tot_xs)*(3.*tau)/(6.-4.*tau+pow(tau,2.));
                                // Upstream angular flux
                                quad_children[i]->psi[j](1,k+1)=\
                                quad_children[i]->psi[j](0,k+1)*(6.-4.*tau)/(6.-4.*tau+pow(tau,2.))+\
                                ((scat_source_lft+ext_source)/tot_xs)*(tau/(6.-4.*tau+pow(tau,2.)))+\
                                ((scat_source_rgt+ext_source)/tot_xs)*(-tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.));
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
        // Directly sweep the directions for the coarsest quadrature
        else {
            dirs_use = quad_selected->get_dirs();
            wgts_use = quad_selected->get_wgts();
            // Go through each direction for the LDFE region
            for (int j=0; j<2; j++) {
                tau = tot_xs*mesh_size/dirs_use(j);
                // Go through each spatial edge
                for (int k=num_edges-2; k>=0; k--) {
                    // Calculate scattering source
                    scat_source_lft=0;
                    scat_source_rgt=0;
                    for (int r=0; r<=scat_order; r++) {
                        leg = legendre(double(r), dirs_use(j));
                        scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k)*leg;
                        scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k)*leg;
                    }
                    // Downstream angular flux
                    quad_selected->psi[j](0,k)=\
                    quad_selected->psi[j](0,k+1)*(6.+2.*tau)/(6.-4.*tau+pow(tau,2.))+\
                    ((scat_source_lft+ext_source)/tot_xs)*(-3.*tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.))-\
                    ((scat_source_rgt+ext_source)/tot_xs)*(3.*tau)/(6.-4.*tau+pow(tau,2.));
                    // Upstream angular flux
                    quad_selected->psi[j](1,k+1)=\
                    quad_selected->psi[j](0,k+1)*(6.-4.*tau)/(6.-4.*tau+pow(tau,2.))+\
                    ((scat_source_lft+ext_source)/tot_xs)*(tau/(6.-4.*tau+pow(tau,2.)))+\
                    ((scat_source_rgt+ext_source)/tot_xs)*(-tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.));
                }
            }
        }
    }
}

double legendre(double order, double ang) {
    
    if (order==0) {
        return 1.;
    }
    else if (order==1) {
        return ang;
    }
    else {
        return ((2.*order-1.)/order)*ang*legendre(order-1.,ang)-((order-1.)/order)*legendre(order-2.,ang);
    }
    
}