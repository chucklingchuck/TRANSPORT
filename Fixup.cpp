//
//  Fixup.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 2/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Fixup.h"

void fixup1 (Region* region_1, Region* region_2, string dir, double tot_source) {
    
    /* Initialize parameters */
    int quad_index;
    int spat_pos_1;
    int spat_pos_2;
    int counter;
    double left_bound;
    double left_dir;
    double right_bound;
    double right_dir;
    double m;
    double b;
    bool reached_bottom;
    vec dir_tmp;
    vec dir_tmp_order;
    uvec sorted_index;
    mat row_tmp_1(1,2);
    mat row_tmp_2(1,2);
    mat interp_table;
    mat interp_table_sorted;
    LDFE_quad* quad_region;
    LDFE_reg* quad_selected;
    vector<LDFE_reg*> quad_ind;
    vector<LDFE_reg*> quad_children;
    vector<LDFE_reg*> quad_children_tmp_ind;
    vector<LDFE_reg*> quad_children_tmp_all;
    vector<mat> psi_tmp;
    
    /* LDFE regions root node index */
    if (dir == "pos") {
        quad_index = 0;
    }
    else {
        quad_index = 1;
    }
    
    /* Spatial edge index */
    if (dir == "pos") {
        spat_pos_1 = region_1->get_num_cells();
        spat_pos_2 = 0;
    }
    else {
        spat_pos_1 = 0;
        spat_pos_2 = region_2->get_num_cells();
    }
    
    /* Find the root node for the region you are mapping from */
    quad_region   = region_1->get_quad();          // This is the quadrature pointer for the region you are mapping from
    quad_ind      = quad_region->get_quad();       // This is the vector of pos and neg LDFE regions for the region you are mapping from
    quad_selected = quad_ind[quad_index];          // This is the positive or negative root LDFE node you have selected
    quad_children = quad_selected->get_children(); // This is the root node's 1st level children
    
    /* Create interpolation table from quadrature you are mapping from */
    // If the quadrature you are mapping from is refined
    if (quad_selected->get_has_children() == true) {
        // Traverse the tree you are mapping from
        counter = 0;
        reached_bottom = false;
        while (reached_bottom == false) {
            // Go through the children nodes in the current level
            for (int i=0; i<quad_children.size(); i++) {
                // If child have children then add them to a temporary children vector
                if (quad_children[i]->get_has_children() == true) {
                    quad_children_tmp_ind = quad_children[i]->get_children();
                    quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                    quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                }
                // If no child then add to interpolation table
                else {
                    psi_tmp = quad_children[i]->psi;
                    dir_tmp = quad_children[i]->get_dirs();
                    row_tmp_1(0,0) = dir_tmp(0);
                    row_tmp_1(0,1) = psi_tmp[0](0,spat_pos_1);
                    row_tmp_2(0,0) = dir_tmp(1);
                    row_tmp_2(0,1) = psi_tmp[1](0,spat_pos_1);
                    interp_table.insert_rows(counter, row_tmp_1);
                    interp_table.insert_rows(counter+1, row_tmp_2);
                    counter = counter+2;
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
    // If the quadrature you are mapping from is the coarsest quadrature
    else if (quad_selected->get_has_children() == false) {
        psi_tmp = quad_selected->psi;
        dir_tmp = quad_selected->get_dirs();
        row_tmp_1(0,0) = dir_tmp(0);
        row_tmp_1(0,1) = psi_tmp[0](0,spat_pos_1);
        row_tmp_2(0,0) = dir_tmp(1);
        row_tmp_2(0,1) = psi_tmp[1](0,spat_pos_1);
        interp_table.insert_rows(0, row_tmp_1);
        interp_table.insert_rows(1, row_tmp_2);
    }
    quad_children.clear();
    
    /* Reorder interpolation table */
    dir_tmp_order.resize(interp_table.n_rows);
    for (int i=0; i<interp_table.n_rows; i++) {
        dir_tmp_order(i) = interp_table(i,0);
    }
    sorted_index = sort_index(dir_tmp_order);
    interp_table_sorted.resize(interp_table.n_rows, 2);
    for (int i=0; i<interp_table.n_rows; i++) {
        interp_table_sorted(i,0) = interp_table(sorted_index(i),0);
        interp_table_sorted(i,1) = interp_table(sorted_index(i),1);
    }
    
    /* Find the root node for the region you are mapping to */
    quad_region   = region_2->get_quad();          // This is the quadrature pointer for the region you are mapping to
    quad_ind      = quad_region->get_quad();       // This is the vector of pos and neg LDFE regions for the region you are mapping to
    quad_selected = quad_ind[quad_index];          // This is the positive or negative root LDFE node you have selected
    quad_children = quad_selected->get_children(); // This is the root node's 1st level children
    
    /* Fixup mapped solution using interpolation table */
    // If the quadrature you are mapping to is refined
    if (quad_selected->get_has_children() == true) {
        // Traverse the tree you are mapping from
        counter = 0;
        reached_bottom = false;
        while (reached_bottom == false) {
            // Go through the children nodes in the current level
            for (int i=0; i<quad_children.size(); i++) {
                // If child have children then add them to a temporary children vector
                if (quad_children[i]->get_has_children() == true) {
                    quad_children_tmp_ind = quad_children[i]->get_children();
                    quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                    quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                }
                // If child has no children then fixup its mapped solution using the interpolation table
                else {
                    psi_tmp = quad_children[i]->psi;
                    dir_tmp = quad_children[i]->get_dirs();
                    // Go through each direction
                    for (int j=0; j<2; j++){
                        // Check if direction is between the leftmost edge and the first direction
                        if (dir_tmp(j) < interp_table_sorted(0,0)) {
                            // For positive directions, leftmost edge is at mu = 0
                            if (dir == "pos") {
                                left_bound = tot_source;
                                left_dir = 0;
                                right_bound = interp_table_sorted(0,1);
                                right_dir = interp_table_sorted(0,0);
                                // Check to see if mapped solution poses an extrema
                                if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                                    (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                                    m = (right_bound-left_bound)/(right_dir-left_dir);
                                    b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                                    quad_children[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                                }
                            }
                            // For negative directions, leftmost edge is at mu = -1
                            else {
                                quad_children[i]->psi[j](0,spat_pos_2) = interp_table_sorted(0,1);
                            }
                        }
                        // Check if direction is between the rightmost edge and the last direction
                        else if (dir_tmp(j) > interp_table_sorted(interp_table_sorted.n_rows-1, 0)) {
                            // For positive directions, rightmost edge is at mu = 1
                            if (dir == "pos") {
                                quad_children[i]->psi[j](0,spat_pos_2) = interp_table_sorted(interp_table_sorted.n_rows-1,1);
                            }
                            // For negative directions, rightmost edge is at mu = 0
                            else {
                                right_bound = tot_source;
                                right_dir = 0;
                                left_bound = interp_table_sorted(interp_table_sorted.n_rows-1,1);
                                left_dir = interp_table_sorted(interp_table_sorted.n_rows-1,0);
                                // Check to see if mapped solution poses an extrema
                                if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                                    (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                                    m = (right_bound-left_bound)/(right_dir-left_dir);
                                    b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                                    quad_children[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                                }
                            }
                        }
                        // Otherwise direction is between two interior directions
                        else {
                            // Go through each row of the sorted interpolation table
                            for (int r=1; r<interp_table_sorted.n_rows; r++) {
                                // Check to see if the correct bounds are found
                                if (dir_tmp(j) < interp_table_sorted(r,0) && dir_tmp(j) > interp_table_sorted(r-1,0)) {
                                    left_dir    = interp_table_sorted(r-1,0);
                                    right_dir   = interp_table_sorted(r,0);
                                    left_bound  = interp_table_sorted(r-1,1);
                                    right_bound = interp_table_sorted(r,1);
                                    // Check to see if mapped solution poses an extrema
                                    if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                                        (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                                        m = (right_bound - left_bound) / (right_dir - left_dir);
                                        b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                                        quad_children[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                                        r=interp_table_sorted.n_rows;
                                    }
                                    // Mapped solution does not need fixup
                                    else {
                                        r=interp_table_sorted.n_rows;
                                    }
                                }
                            }
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
    // If the quadrature you are mapping to is the coarsest quadrature
    else if (quad_selected->get_has_children() == false) {
        psi_tmp = quad_selected->psi;
        dir_tmp = quad_selected->get_dirs();
        // Look up each direction
        for (int j=0; j<2; j++){
            // Check if direction is between the leftmost edge and the first direction
            if (dir_tmp(j) < interp_table_sorted(0,0)) {
                // For positive directions, leftmost edge is at mu=0
                if (dir == "pos") {
                    left_bound = tot_source;
                    left_dir = 0;
                    right_bound = interp_table_sorted(0,1);
                    right_dir = interp_table_sorted(0,0);
                    // Check to see if mapped solution poses an extrema
                    if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                        (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                        m = (right_bound-left_bound)/(right_dir-left_dir);
                        b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                        quad_selected->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                    }
                }
                // For negative directions, the leftmost edge is at mu=-1
                else {
                    quad_selected->psi[j](0,spat_pos_2) = interp_table_sorted(0,1);
                }
            }
            // Check if direction is between the rightmost edge and the last direction
            else if (dir_tmp(j) > interp_table_sorted(interp_table_sorted.n_rows-1, 0)) {
                // For positive directions, rightmost edge is at mu=1
                if (dir == "pos") {
                    quad_selected->psi[j](0,spat_pos_2) = interp_table_sorted(interp_table_sorted.n_rows-1,1);
                }
                // For negative directions, rightmost edge is at mu=0
                else {
                    right_bound = tot_source;
                    right_dir = 0;
                    left_bound = interp_table_sorted(interp_table_sorted.n_rows-1,1);
                    left_dir = interp_table_sorted(interp_table_sorted.n_rows-1,0);
                    // Check to see if mapped solution poses an extrema
                    if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                        (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                        m = (right_bound - left_bound) / (right_dir - left_dir);
                        b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                        quad_selected->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                    }
                }
            }
            // Otherwise direction is between two interior directions
            else {
                // Go through each row of the sorted interpolation table
                for (int r=1; r<interp_table_sorted.n_rows; r++) {
                    // Check to see if the correct bounds are found
                    if (dir_tmp(j) < interp_table_sorted(r,0) && dir_tmp(j) > interp_table_sorted(r-1,0)) {
                        right_bound = interp_table_sorted(r,1);
                        right_dir = interp_table_sorted(r,0);;
                        left_bound = interp_table_sorted(r-1,1);
                        left_dir = interp_table_sorted(r-1,0);
                        // Check to see if mapped solution poses an extrema
                        if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
                            (psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
                            m = (right_bound-left_bound)/(right_dir-left_dir);
                            b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
                            quad_selected->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
                        }
                        // Mapped solution does not need fixup
                        else {
                            r=interp_table_sorted.n_rows;
                        }
                    }
                }
            }
        }
    }
    quad_children.clear();
}