//
//  Mapping.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Mapping.h"

void mapping1 (Region* region_1, Region* region_2, string dir) {
    
    /* Initialize parameters */
    int quad_index;
    int spat_pos_1;
    int spat_pos_2;
    bool reached_bottom;
    bool reached_bottom_inner;
    LDFE_quad* quad_region_1;
    LDFE_quad* quad_region_2;
    LDFE_reg* quad_selected_1;
    LDFE_reg* quad_selected_2;
    vector<LDFE_reg*> quad_ind_1;
    vector<LDFE_reg*> quad_ind_2;
    vector<LDFE_reg*> quad_children_1;
    vector<LDFE_reg*> quad_children_2;
    vector<LDFE_reg*> quad_children_tmp_ind_1;
    vector<LDFE_reg*> quad_children_tmp_ind_2;
    vector<LDFE_reg*> quad_children_tmp_all_1;
    vector<LDFE_reg*> quad_children_tmp_all_2;
    vector<LDFE_reg*> quad_inner_children;
    vector<mat> map_values;
    vector<LDFE_reg*> quad_inner_children_tmp_ind;
    vector<LDFE_reg*> quad_inner_children_tmp_all;
    mat map_tmp(2, 2);
    int level;
    int col_actual;
    double ave_left;
    double ave_right;
    bool found;
    bool at_beginning;
    int col_loc;
    
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
    quad_region_1   = region_1->get_quad();            // Quadrature pointer for the region you are mapping from
    quad_ind_1      = quad_region_1->get_quad();       // Vector of pos and neg LDFE regions for the region you are mapping from
    quad_selected_1 = quad_ind_1[quad_index];          // Positive or negative root LDFE node you have selected
    quad_children_1 = quad_selected_1->get_children(); // Root node's 1st level children
    
    /* Find the root node for the region you are mapping to */
    quad_region_2   = region_2->get_quad();            // Quadrature pointer for the region you are mapping to
    quad_ind_2      = quad_region_2->get_quad();       // Vector of pos and neg LDFE regions for the region you are mapping to
    quad_selected_2 = quad_ind_2[quad_index];          // Positive or negative root LDFE node you have selected
    quad_children_2 = quad_selected_2->get_children(); // Root node's 1st level children
    
    /* Case 1: Both quadratures have children */
    if (quad_selected_1->get_has_children() == true && quad_selected_2->get_has_children() == true) {
        // Traverse the tree you are mapping from while comparing to the tree you are mapping to
        reached_bottom = false;
        while (reached_bottom == false) {
            // Go through the children nodes in the current level
            for (int i=0; i<quad_children_1.size(); i++) {
                // If both children have children then add them to a temporary children vector
                if (quad_children_1[i]->get_has_children() == true && quad_children_2[i]->get_has_children() == true) {
                    quad_children_tmp_ind_1 = quad_children_1[i]->get_children();
                    quad_children_tmp_all_1.push_back(quad_children_tmp_ind_1[0]);
                    quad_children_tmp_all_1.push_back(quad_children_tmp_ind_1[1]);
                    quad_children_tmp_ind_2 = quad_children_2[i]->get_children();
                    quad_children_tmp_all_2.push_back(quad_children_tmp_ind_2[0]);
                    quad_children_tmp_all_2.push_back(quad_children_tmp_ind_2[1]);
                }
                // Fine-to-coarse mapping
                else if (quad_children_1[i]->get_has_children() == true && quad_children_2[i]->get_has_children() == false) {
                    // Set the child node as the inner root node
                    quad_inner_children = quad_children_1[i]->get_children(); // Inner root node's first level children
                    // Traverse to the bottom of the inner root node
                    level = 1; // Set first level as one. The final mapped solution should be at the zeroth level
                    reached_bottom_inner = false;
                    while (reached_bottom_inner == false) {
                        for (int k=0; k<quad_inner_children.size(); k++) { // Go through the inner children nodes in the current level
                            // If the current child has children add them to a temporary children vector
                            if (quad_inner_children[k]->get_has_children() == true) {
                                quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                            }
                            // Store the psi, level and column position of nodes without children
                            else {
                                for (int r=0; r<2; r++) {
                                    map_tmp(1,r) = quad_inner_children[k]->psi[r](0,spat_pos_1);
                                }
                                map_tmp(0,0) = level;
                                map_tmp(0,1) = k;
                                map_values.push_back(map_tmp);
                            }
                        }
                        // Check if the bottom of the inner tree has been reached
                        if (quad_inner_children_tmp_all.size() == 0) {
                            reached_bottom_inner = true;
                        }
                        else {
                            level++;
                        }
                        // Set the inner children vector as the temporary children vetor for the next level calculation
                        quad_inner_children = quad_inner_children_tmp_all;
                        quad_inner_children_tmp_all.clear();
                        quad_inner_children_tmp_ind.clear();
                    }
                    // Mapped stored psi back to child node
                    // Go through all the levels from the bottom up
                    for (int k=level; k>0; k--) {
                        // Go through the entire map_values vector
                        for (int r=0; r<map_values.size(); r++) {
                            // Find the first element in map_values at the current level
                            if (int(map_values[r](0,0))==k) {
                                // Average the angular fluxes in the current and adjacent LDFE ranges
                                ave_left = (map_values[r](1,0)+map_values[r](1,1))/2;
                                ave_right = (map_values[r+1](1,0)+map_values[r+1](1,1))/2;
                                // Delete the two elements that were averaged together
                                map_values.erase(map_values.begin()+r, map_values.begin()+r+2);
                                // Find col_actual which is the position in map_values to insert the new map_tmp
                                // Find col_loc which is the position relative to other elements in the same level
                                found = false;
                                at_beginning = true;
                                col_loc = 0; // Start searching for the leftmost element in the same level
                                // Go through the entire map_values vector
                                for (int t=0; t<map_values.size(); t++) {
                                    // Find the elements at the level of the new map_tmp
                                    if (map_values[t](0,0)==k-1) {
                                        // If such an element is found then col_actual cannot be at the beginning of map_values
                                        at_beginning = false;
                                        // Check if the element is at the current col_loc
                                        if (map_values[t](0,1)==col_loc) {
                                            col_loc++;        // Incrememnt col_loc
                                            col_actual = t+1; // Set col_actual as the adjacent position
                                            found = true;
                                        }
                                        // Check if the element should be inserted between two elements at the level of the new map_tmp
                                        else if (map_values[t](0,1)>col_loc) {
                                            col_actual = t;             // Set col_actual as the current position
                                            t = int(map_values.size()); // Exit loop
                                            found = true;
                                        }
                                    }
                                }
                                // The new map_tmp should be at the beginning of the map_values vector
                                if (at_beginning == true) {
                                    col_actual = 0;
                                }
                                // The new map_tmp should be at the end of the map_values vector
                                else if (found == false) {
                                    col_actual = int(map_values.size()) ;
                                }
                                // Create the new map_tmp
                                map_tmp(0,0) = k-1;
                                map_tmp(0,1) = col_loc;
                                map_tmp(1,0) = ave_left;
                                map_tmp(1,1) = ave_right;
                                // Insert the new map_tmp into map_values
                                map_values.insert(map_values.begin()+col_actual, map_tmp);
                                r=0; // Reset r
                            }
                        }
                    }
                    map_values.clear();
                    // Pass the mapped values to the downstream LDFE region
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0,spat_pos_2) = map_tmp(1,k);
                    }
                    quad_inner_children.clear();
                }
                // Coarse-to-fine mapping
                else if (quad_children_1[i]->get_has_children() == false && quad_children_2[i]->get_has_children() == true) {
                    // Pass solution to quadrature you are mapping to
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0, spat_pos_2) = quad_children_1[i]->psi[k](0, spat_pos_1);
                    }
                    // Set the child node as the root node
                    quad_inner_children = quad_children_2[i]->get_children(); // Root node's first level of children
                    // Map to the first level children
                    quad_inner_children[0]->psi.resize(2);
                    quad_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[0]->psi[0](0, spat_pos_2) = quad_children_2[i]->psi[0](0, spat_pos_2);
                    quad_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[0]->psi[1](0, spat_pos_2) = quad_children_2[i]->psi[0](0, spat_pos_2);
                    quad_inner_children[1]->psi.resize(2);
                    quad_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[1]->psi[0](0, spat_pos_2) = quad_children_2[i]->psi[1](0, spat_pos_2);
                    quad_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[1]->psi[1](0, spat_pos_2) = quad_children_2[i]->psi[1](0, spat_pos_2);
                    // Traverse to the bottom of the root node
                    reached_bottom_inner = false;
                    while (reached_bottom_inner == false) {
                        for (int k=0; k<quad_inner_children.size(); k++) {
                            // If the current child has children add them to a temporary children vector
                            if (quad_inner_children[k]->get_has_children() == true) {
                                quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                                // Add the parent flux to children nodes
                                vector<LDFE_reg*> quad_inner_inner_children = quad_inner_children[k]->get_children();
                                quad_inner_inner_children[0]->psi.resize(2);
                                quad_inner_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[0]->psi[0](0, spat_pos_2) = quad_inner_children[k]->psi[0](0, spat_pos_2);
                                quad_inner_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[0]->psi[1](0, spat_pos_2) = quad_inner_children[k]->psi[0](0, spat_pos_2);
                                quad_inner_inner_children[1]->psi.resize(2);
                                quad_inner_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[1]->psi[0](0, spat_pos_2) = quad_inner_children[k]->psi[1](0, spat_pos_2);
                                quad_inner_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[1]->psi[1](0, spat_pos_2) = quad_inner_children[k]->psi[1](0, spat_pos_2);
                            }
                        }
                        // Check if the bottom of the inner tree has been reached
                        if (quad_inner_children_tmp_all.size() == 0) {
                            reached_bottom_inner = true;
                        }
                        // Set the inner children vector as the temporary children vetor for the next calculation
                        quad_inner_children = quad_inner_children_tmp_all;
                        quad_inner_children_tmp_all.clear();
                        quad_inner_children_tmp_ind.clear();
                    }
                    quad_inner_children.clear();
                }
                // One-to-one mapping if both children do not have children
                else {
                    /* Pass the solution to the quadrature you are mapping to. */
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0, spat_pos_2) = quad_children_1[i]->psi[k](0, spat_pos_1);
                    }
                }
            }
            // If there are no children in the temporary children vector the bottom is reached
            if (quad_children_tmp_all_1.size() == 0 && quad_children_tmp_all_2.size() == 0) {
                reached_bottom = true;
            }
            // Set the children vector as the temporary children vector for the next level calculation
            quad_children_1 = quad_children_tmp_all_1;
            quad_children_2 = quad_children_tmp_all_2;
            quad_children_tmp_all_1.clear();
            quad_children_tmp_all_2.clear();
            quad_children_tmp_ind_1.clear();
            quad_children_tmp_ind_2.clear();
        }
    }
    
    /* Case 2: Mapping from a refined quadrature to a base quadrature */
    else if (quad_selected_1->get_has_children() == true && quad_selected_2->get_has_children() == false) {
        // Traverse to the bottom of the inner root node
        level = 1; // Set first level as one. The final mapped solution should be at the zeroth level
        quad_inner_children  = quad_children_1;
        reached_bottom = false;
        while (reached_bottom == false) {
            for (int k=0; k<quad_inner_children.size(); k++) { // Go through the inner children nodes in the current level
                // If the current child has children add them to a temporary children vector
                if (quad_inner_children[k]->get_has_children() == true) {
                    quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                }
                // Store the psi, level and column position of nodes without children
                else {
                    for (int r=0; r<2; r++) {
                        map_tmp(1,r) = quad_inner_children[k]->psi[r](0,spat_pos_1);
                    }
                    map_tmp(0,0) = level;
                    map_tmp(0,1) = k;
                    map_values.push_back(map_tmp);
                }
            }
            // Check if the bottom of the inner tree has been reached
            if (quad_inner_children_tmp_all.size() == 0) {
                reached_bottom = true;
            }
            else {
                level++;
            }
            // Set the inner children vector as the temporary children vector for the next level calculation
            quad_inner_children = quad_inner_children_tmp_all;
            quad_inner_children_tmp_all.clear();
            quad_inner_children_tmp_ind.clear();
        }
        // Mapped stored psi back to child node
        // Go through all the levels from the bottom up
        for (int k=level; k>0; k--) {
            // Go through the entire map_values vector
            for (int r=0; r<map_values.size(); r++) {
                // Find the first element in map_values at the current level
                if (int(map_values[r](0,0))==k) {
                    // Average the angular fluxes in the current and adjacent LDFE ranges
                    ave_left = (map_values[r](1,0)+map_values[r](1,1))/2;
                    ave_right = (map_values[r+1](1,0)+map_values[r+1](1,1))/2;
                    // Delete the two elements that were averaged together
                    map_values.erase(map_values.begin()+r, map_values.begin()+r+2);
                    // Find col_actual which is the position in map_values to insert the new map_tmp
                    // Find col_loc which is the position relative to other elements in the same level
                    found = false;
                    at_beginning = true;
                    col_loc = 0; // Start searching for the leftmost element in the same level
                    // Go through the entire map_values vector
                    for (int t=0; t<map_values.size(); t++) {
                        // Find the elements at the level of the new map_tmp
                        if (map_values[t](0,0)==k-1) {
                            // If such an element is found then col_actual cannot be at the beginning of map_values
                            at_beginning = false;
                            // Check if the element is at the current col_loc
                            if (map_values[t](0,1)==col_loc) {
                                col_loc++;        // Increment col_loc
                                col_actual = t+1; // Set col_actual as the adjacent position
                                found = true;
                            }
                            // Check if the element should be inserted between two elements at the level of the new map_tmp
                            else if (map_values[t](0,1)>col_loc) {
                                col_actual = t;             // Set col_actual as the current position
                                t = int(map_values.size()); // Exit loop
                                found = true;
                            }
                        }
                    }
                    // The new map_tmp should be at the beginning of the map_values vector
                    if (at_beginning == true) {
                        col_actual = 0;
                    }
                    // The new map_tmp should be at the end of the map_values vector
                    else if (found == false) {
                        col_actual = int(map_values.size()) ;
                    }
                    // Create the new map_tmp
                    map_tmp(0,0) = k-1;
                    map_tmp(0,1) = col_loc;
                    map_tmp(1,0) = ave_left;
                    map_tmp(1,1) = ave_right;
                    // Insert the new map_tmp into map_values
                    map_values.insert(map_values.begin()+col_actual, map_tmp);
                    r=0; // Reset r
                }
            }
        }
        map_values.clear();
        // Pass the mapped values to the downstream LDFE region
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0,spat_pos_2) = map_tmp(1,k);
        }
        quad_inner_children.clear();
    }
    
    /* Case 3: Mapping from a base quadrature to a refined quadrature */
    else if (quad_selected_1->get_has_children() == false && quad_selected_2->get_has_children() == true) {
        // Pass solution to quadrature you are mapping to
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0, spat_pos_2) = quad_selected_1->psi[k](0, spat_pos_1);
        }
        // Set the child node as the root node
        quad_inner_children = quad_selected_2->get_children(); // Root node's first level of children
        // Map to the first level children
        quad_inner_children[0]->psi.resize(2);
        quad_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[0]->psi[0](0, spat_pos_2) = quad_selected_2->psi[0](0, spat_pos_2);
        quad_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[0]->psi[1](0, spat_pos_2) = quad_selected_2->psi[0](0, spat_pos_2);
        quad_inner_children[1]->psi.resize(2);
        quad_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[1]->psi[0](0, spat_pos_2) = quad_selected_2->psi[1](0, spat_pos_2);
        quad_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[1]->psi[1](0, spat_pos_2) = quad_selected_2->psi[1](0, spat_pos_2);
        // Traverse to the bottom of the root node
        reached_bottom_inner = false;
        while (reached_bottom_inner == false) {
            for (int k=0; k<quad_inner_children.size(); k++) {
                // If the current child has children add them to a temporary children vector
                if (quad_inner_children[k]->get_has_children() == true) {
                    quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                    // Add the parent flux to children nodes
                    vector<LDFE_reg*> quad_inner_inner_children = quad_inner_children[k]->get_children();
                    quad_inner_inner_children[0]->psi.resize(2);
                    quad_inner_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[0]->psi[0](0, spat_pos_2) = quad_inner_children[k]->psi[0](0, spat_pos_2);
                    quad_inner_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[0]->psi[1](0, spat_pos_2) = quad_inner_children[k]->psi[0](0, spat_pos_2);
                    quad_inner_inner_children[1]->psi.resize(2);
                    quad_inner_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[1]->psi[0](0, spat_pos_2) = quad_inner_children[k]->psi[1](0, spat_pos_2);
                    quad_inner_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[1]->psi[1](0, spat_pos_2) = quad_inner_children[k]->psi[1](0, spat_pos_2);
                }
            }
            // Check if the bottom of the inner tree has been reached
            if (quad_inner_children_tmp_all.size() == 0) {
                reached_bottom_inner = true;
            }
            // Set the inner children vector as the temporary children vetor for the next calculation
            quad_inner_children = quad_inner_children_tmp_all;
            quad_inner_children_tmp_all.clear();
            quad_inner_children_tmp_ind.clear();
        }
        quad_inner_children.clear();
    }
    
    /* Mapping between two base quadratures */
    else {
        // Pass the solution to the quadrature you are mapping to
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0, spat_pos_2) = quad_selected_1->psi[k](0, spat_pos_1);
        }
    }
}

void mapping2 (Region* region_1, Region* region_2, string dir) {
    
    /* Initialize parameters */
    int quad_index;
    int spat_pos_1;
    int spat_pos_2;
    int counter;
    int level;
    int col_actual;
    int col_loc;
    bool reached_bottom;
    bool reached_bottom_inner;
    bool found;
    bool at_beginning;
    LDFE_quad* quad_region_1;
    LDFE_quad* quad_region_2;
    LDFE_reg* quad_selected_1;
    LDFE_reg* quad_selected_2;
    vector<LDFE_reg*> quad_ind_1;
    vector<LDFE_reg*> quad_ind_2;
    vector<LDFE_reg*> quad_children_1;
    vector<LDFE_reg*> quad_children_2;
    vector<LDFE_reg*> quad_children_tmp_ind_1;
    vector<LDFE_reg*> quad_children_tmp_ind_2;
    vector<LDFE_reg*> quad_children_tmp_all_1;
    vector<LDFE_reg*> quad_children_tmp_all_2;
    vector<LDFE_reg*> quad_children_tmp_ind_map;
    vector<LDFE_reg*> quad_children_tmp_all_map;
    vector<LDFE_reg*> quad_children_map;
    vector<LDFE_reg*> quad_inner_children;
    vector<LDFE_reg*> quad_inner_children_tmp_ind;
    vector<LDFE_reg*> quad_inner_children_tmp_all;
    vector<vector<mat> > basis_tmp;
    vector<vector<vec> > wgt_tmp;
    vector<vector<vec> > dir_tmp;
    vector<mat> map_values;
    vector<mat> sol_tmp_coarse;
    vec b_vec(2);
    vec dirs_tmp_ind;
    vec dir_wgts_ind;
    vec coarse_sol;
    vec wgts_tmp_coarse;
    vec dirs_tmp_coarse;
    vec dirs_tmp_fine_0;
    vec dirs_tmp_fine_1;
    vec wgts_tmp_fine_0;
    vec wgts_tmp_fine_1;
    vec coarse_tilde;
    vec fine_sol_0(2);
    vec fine_sol_1(2);
    vec dirs_check;
    mat basis_ind;
    mat a_mat(2,2);
    mat map_tmp(6, 2);
    mat basis_tmp_coarse;
    mat basis_tmp_fine_0;
    mat basis_tmp_fine_1;
    
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
    quad_region_1   = region_1->get_quad();            // Quadrature pointer for the region you are mapping from
    quad_ind_1      = quad_region_1->get_quad();       // Vector of pos and neg LDFE regions for the region you are mapping from
    quad_selected_1 = quad_ind_1[quad_index];          // Positive or negative root LDFE node you have selected
    quad_children_1 = quad_selected_1->get_children(); // Root node's 1st level children
    
    /* Find the root node for the region you are mapping to */
    quad_region_2   = region_2->get_quad();            // Quadrature pointer for the region you are mapping to
    quad_ind_2      = quad_region_2->get_quad();       // Vector of pos and neg LDFE regions for the region you are mapping to
    quad_selected_2 = quad_ind_2[quad_index];          // Positive or negative root LDFE node you have selected
    quad_children_2 = quad_selected_2->get_children(); // Root node's 1st level children
    
    /* Case 1: Both quadratures have children */
    if (quad_selected_1->get_has_children() == true && quad_selected_2->get_has_children() == true) {
        // Traverse the tree you are mapping from while comparing to the tree you are mapping to
        reached_bottom = false;
        while (reached_bottom == false) {
            // Go through the children nodes in the current level
            for (int i=0; i<quad_children_1.size(); i++) {
                // If both childs have children then add them to a temporary children vector
                if (quad_children_1[i]->get_has_children() == true && quad_children_2[i]->get_has_children() == true) {
                    quad_children_tmp_ind_1 = quad_children_1[i]->get_children();
                    quad_children_tmp_all_1.push_back(quad_children_tmp_ind_1[0]);
                    quad_children_tmp_all_1.push_back(quad_children_tmp_ind_1[1]);
                    quad_children_tmp_ind_2 = quad_children_2[i]->get_children();
                    quad_children_tmp_all_2.push_back(quad_children_tmp_ind_2[0]);
                    quad_children_tmp_all_2.push_back(quad_children_tmp_ind_2[1]);
                }
                // Fine-to-coarse mapping
                else if (quad_children_1[i]->get_has_children() == true && quad_children_2[i]->get_has_children() == false) {
                    // Set the child node as the inner root node
                    quad_inner_children = quad_children_1[i]->get_children(); // Inner root node's first level children
                    // Add the basis and quadrature of the root node
                    basis_tmp.push_back(vector<mat>());
                    wgt_tmp.push_back(vector<vec>());
                    dir_tmp.push_back(vector<vec>());
                    basis_tmp[0].push_back(quad_children_1[i]->get_basis());
                    wgt_tmp[0].push_back(quad_children_1[i]->get_wgts());
                    dir_tmp[0].push_back(quad_children_1[i]->get_dirs());
                    // Add space to store information for the next level of children
                    basis_tmp.push_back(vector<mat>());
                    wgt_tmp.push_back(vector<vec>());
                    dir_tmp.push_back(vector<vec>());
                    // Traverse to the bottom of the inner root node
                    level = 1;
                    reached_bottom_inner = false;
                    while (reached_bottom_inner == false) {
                        // Go through the inner children nodes in the current level
                        for (int k=0; k<quad_inner_children.size(); k++) {
                            // If the current child has children add them to a temporary children vector
                            if (quad_inner_children[k]->get_has_children() == true) {
                                quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                                // Store the basis and quadrature for mapping
                                basis_tmp[level].push_back(quad_inner_children[k]->get_basis());
                                wgt_tmp[level].push_back(quad_inner_children[k]->get_wgts());
                                dir_tmp[level].push_back(quad_inner_children[k]->get_dirs());
                            }
                            // Store the psi, level, column position, basis and quadrature of nodes without children
                            else {
                                // Store level and column position
                                map_tmp(0,0) = level;
                                map_tmp(0,1) = k;
                                // Store psi
                                for (int r=0; r<2; r++) {
                                    map_tmp(1,r) = quad_inner_children[k]->psi[r](0,spat_pos_1);
                                }
                                // Store directions
                                dirs_tmp_ind = quad_inner_children[k]->get_dirs();
                                map_tmp(2,0) = dirs_tmp_ind(0);
                                map_tmp(2,1) = dirs_tmp_ind(1);
                                // Store weights
                                dir_wgts_ind = quad_inner_children[k]->get_wgts();
                                map_tmp(3,0) = dir_wgts_ind(0);
                                map_tmp(3,1) = dir_wgts_ind(1);
                                // Store basis
                                basis_ind = quad_inner_children[k]->get_basis();
                                map_tmp(4,0) = basis_ind(0,0);
                                map_tmp(4,1) = basis_ind(0,1);
                                map_tmp(5,0) = basis_ind(1,0);
                                map_tmp(5,1) = basis_ind(1,1);
                                map_values.push_back(map_tmp);
                            }
                        }
                        // Check if the bottom of the inner tree has been reached
                        if (quad_inner_children_tmp_all.size() == 0) {
                            reached_bottom_inner = true;
                        }
                        else {
                            level++;
                            basis_tmp.push_back(vector<mat>());
                            wgt_tmp.push_back(vector<vec>());
                            dir_tmp.push_back(vector<vec>());
                        }
                        // Set the inner children vector as the temporary children vetor for the next level calculation
                        quad_inner_children = quad_inner_children_tmp_all;
                        quad_inner_children_tmp_all.clear();
                        quad_inner_children_tmp_ind.clear();
                    }
                    // Mapped stored psi back to child node
                    // Go through all the levels from the bottom up
                    for (int k=level; k>0; k--) {
                        // Track coarse nodes in current level
                        counter = 0;
                        // Go through the map_values vector
                        for (int r=0; r<map_values.size(); r++) {
                            // Find the first element in map_values at the current level
                            if (int(map_values[r](0,0))==k) {
                                // Create the b-vector
                                b_vec[0] = map_values[r](3,0)*map_values[r](1,0) + \
                                map_values[r](3,1)*map_values[r](1,1) + \
                                map_values[r+1](3,0)*map_values[r+1](1,0) + \
                                map_values[r+1](3,1)*map_values[r+1](1,1);
                                b_vec[1] = map_values[r](2,0)*map_values[r](3,0)*map_values[r](1,0) + \
                                map_values[r](2,1)*map_values[r](3,1)*map_values[r](1,1) + \
                                map_values[r+1](2,0)*map_values[r+1](3,0)*map_values[r+1](1,0) + \
                                map_values[r+1](2,1)*map_values[r+1](3,1)*map_values[r+1](1,1);
                                // Create the a-matrix
                                a_mat(0,0) = (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](0))*wgt_tmp[k-1][counter](0) + \
                                (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](1))*wgt_tmp[k-1][counter](1);
                                a_mat(0,1) = (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](0))*wgt_tmp[k-1][counter](0) + \
                                (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](1))*wgt_tmp[k-1][counter](1);
                                a_mat(1,0) = (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](0))*dir_tmp[k-1][counter](0)*wgt_tmp[k-1][counter](0) + \
                                (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](1))*dir_tmp[k-1][counter](1)*wgt_tmp[k-1][counter](1);
                                a_mat(1,1) = (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](0))*dir_tmp[k-1][counter](0)*wgt_tmp[k-1][counter](0) + \
                                (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](1))*dir_tmp[k-1][counter](1)*wgt_tmp[k-1][counter](1);
                                coarse_sol = solve(a_mat, b_vec);
                                // Delete the two elements that were averaged together
                                map_values.erase(map_values.begin()+r, map_values.begin()+r+2);
                                // Find col_actual which is the position in map_values to insert the new map_tmp
                                // Find col_loc which is the position relative to other elements in the same level
                                found = false;
                                at_beginning = true;
                                // Start searching for the leftmost element in the same level
                                col_loc = 0;
                                // Go through the entire map_values vector
                                for (int t=0; t<map_values.size(); t++) {
                                    // Find the elements at the level of the new map_tmp
                                    if (map_values[t](0,0)==k-1) {
                                        // If such an element is found then col_actual cannot be at the beginning of map_values
                                        at_beginning = false;
                                        // Check if the element is at the current col_loc
                                        if (map_values[t](0,1)==col_loc) {
                                            col_loc++;        // Incrememnt col_loc
                                            col_actual = t+1; // Set col_actual as the adjacent position
                                            found = true;
                                        }
                                        // Check if the element should be inserted between two elements at the level of the new map_tmp
                                        else if (map_values[t](0,1)>col_loc) {
                                            col_actual = t;             // Set col_actual as the current position
                                            t = int(map_values.size()); // Exit loop
                                            found = true;
                                        }
                                    }
                                }
                                // The new map_tmp should be at the beginning of the map_values vector
                                if (at_beginning == true) {
                                    col_actual = 0;
                                }
                                // The new map_tmp should be at the end of the map_values vector
                                else if (found == false) {
                                    col_actual = int(map_values.size()) ;
                                }
                                // Create the new map_tmp
                                map_tmp(0,0) = k-1;
                                map_tmp(0,1) = col_loc;
                                map_tmp(1,0) = coarse_sol(0);
                                map_tmp(1,1) = coarse_sol(1);
                                map_tmp(2,0) = dir_tmp[k-1][counter](0);
                                map_tmp(2,1) = dir_tmp[k-1][counter](1);
                                map_tmp(3,0) = wgt_tmp[k-1][counter](0);
                                map_tmp(3,1) = wgt_tmp[k-1][counter](1);
                                map_tmp(4,0) = basis_tmp[k-1][counter](0,0);
                                map_tmp(4,1) = basis_tmp[k-1][counter](0,1);
                                map_tmp(5,0) = basis_tmp[k-1][counter](1,0);
                                map_tmp(5,1) = basis_tmp[k-1][counter](1,1);
                                // Insert the new map_tmp into map_values
                                map_values.insert(map_values.begin()+col_actual, map_tmp);
                                // Reset r
                                r=0;
                                // Increment counter to go to next coarse node in current level
                                counter = counter + 1;
                            }
                        }
                    }
                    map_values.clear();
                    basis_tmp.clear();
                    wgt_tmp.clear();
                    dir_tmp.clear();
                    // Pass the mapped values to the downstream LDFE region
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0,spat_pos_2) = map_tmp(1,k);
                    }
                    quad_inner_children.clear();
                }
                // Coarse-to-fine mapping
                else if (quad_children_1[i]->get_has_children() == false && quad_children_2[i]->get_has_children() == true) {
                    // Pass solution to quadrature you are mapping to
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0, spat_pos_2) = quad_children_1[i]->psi[k](0, spat_pos_1);
                    }
                    // Set the child node as the root node
                    quad_inner_children = quad_children_2[i]->get_children(); // Root node's first level of children
                    // Map to the first level children
                    basis_tmp_coarse = quad_children_2[i]->get_basis();
                    wgts_tmp_coarse = quad_children_2[i]->get_wgts();
                    dirs_tmp_coarse = quad_children_2[i]->get_dirs();
                    sol_tmp_coarse = quad_children_2[i]->psi;
                    dirs_tmp_fine_0 = quad_inner_children[0]->get_dirs();
                    dirs_tmp_fine_1 = quad_inner_children[1]->get_dirs();
                    wgts_tmp_fine_0 = quad_inner_children[0]->get_wgts();
                    wgts_tmp_fine_1 = quad_inner_children[1]->get_wgts();
                    basis_tmp_fine_0 = quad_inner_children[0]->get_basis();
                    basis_tmp_fine_1 = quad_inner_children[1]->get_basis();
                    // Assemble A-matrix
                    a_mat(0,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1));
                    a_mat(0,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                    a_mat(1,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                    a_mat(1,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                    // Assemble b-vector
                    b_vec(0) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)+\
                    wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2);
                    b_vec(1) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)*dirs_tmp_coarse(0)+\
                    wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2)*dirs_tmp_coarse(1);
                    // Solve for coarse tilde
                    coarse_tilde = solve(a_mat, b_vec);
                    // Solve for fine psi
                    fine_sol_0(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0));
                    fine_sol_0(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1));
                    fine_sol_1(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0));
                    fine_sol_1(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                    // Store fine psi into next level children
                    quad_inner_children[0]->psi.resize(2);
                    quad_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[0]->psi[0](0, spat_pos_2) = fine_sol_0(0);
                    quad_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[0]->psi[1](0, spat_pos_2) = fine_sol_0(1);
                    quad_inner_children[1]->psi.resize(2);
                    quad_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[1]->psi[0](0, spat_pos_2) = fine_sol_1(0);
                    quad_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_children[1]->psi[1](0, spat_pos_2) = fine_sol_1(1);
                    // Traverse to the bottom of the root node
                    reached_bottom_inner = false;
                    while (reached_bottom_inner == false) {
                        for (int k=0; k<quad_inner_children.size(); k++) {
                            // If the current child has children add them to a temporary children vector
                            if (quad_inner_children[k]->get_has_children() == true) {
                                quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                                quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                                // Add the parent flux to children nodes
                                vector<LDFE_reg*> quad_inner_inner_children = quad_inner_children[k]->get_children();
                                // Map to the first level children
                                basis_tmp_coarse = quad_inner_children[k]->get_basis();
                                wgts_tmp_coarse = quad_inner_children[k]->get_wgts();
                                dirs_tmp_coarse = quad_inner_children[k]->get_dirs();
                                sol_tmp_coarse = quad_inner_children[k]->psi;
                                dirs_tmp_fine_0 = quad_inner_inner_children[0]->get_dirs();
                                dirs_tmp_fine_1 = quad_inner_inner_children[1]->get_dirs();
                                wgts_tmp_fine_0 = quad_inner_inner_children[0]->get_wgts();
                                wgts_tmp_fine_1 = quad_inner_inner_children[1]->get_wgts();
                                basis_tmp_fine_0 = quad_inner_inner_children[0]->get_basis();
                                basis_tmp_fine_1 = quad_inner_inner_children[1]->get_basis();
                                // Assemble A-matrix
                                a_mat(0,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                                wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                                wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                                wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1));
                                a_mat(0,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))+\
                                wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))+\
                                wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))+\
                                wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                                a_mat(1,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                                wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                                wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                                wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                                a_mat(1,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                                wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                                wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                                wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                                // Assemble b-vector
                                b_vec(0) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)+\
                                wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2);
                                b_vec(1) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)*dirs_tmp_coarse(0)+\
                                wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2)*dirs_tmp_coarse(1);
                                // Solve for coarse tilde
                                coarse_tilde = solve(a_mat, b_vec);
                                // Solve for fine psi
                                fine_sol_0(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                                coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0));
                                fine_sol_0(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                                coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1));
                                fine_sol_1(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                                coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0));
                                fine_sol_1(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))+\
                                coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                                // Store fine psi into next level children
                                quad_inner_inner_children[0]->psi.resize(2);
                                quad_inner_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[0]->psi[0](0, spat_pos_2) = fine_sol_0(0);
                                quad_inner_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[0]->psi[1](0, spat_pos_2) = fine_sol_0(1);
                                quad_inner_inner_children[1]->psi.resize(2);
                                quad_inner_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[1]->psi[0](0, spat_pos_2) = fine_sol_1(0);
                                quad_inner_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                                quad_inner_inner_children[1]->psi[1](0, spat_pos_2) = fine_sol_1(1);
                            }
                        }
                        // Check if the bottom of the inner tree has been reached
                        if (quad_inner_children_tmp_all.size() == 0) {
                            reached_bottom_inner = true;
                        }
                        // Set the inner children vector as the temporary children vetor for the next calculation
                        quad_inner_children = quad_inner_children_tmp_all;
                        quad_inner_children_tmp_all.clear();
                        quad_inner_children_tmp_ind.clear();
                    }
                    quad_inner_children.clear();
                }
                // One-to-one mapping if both children do not have children
                else {
                    /* Pass the solution to the quadrature you are mapping to. */
                    quad_children_2[i]->psi.resize(2);
                    for (int k=0; k<2; k++) {
                        quad_children_2[i]->psi[k].resize(2, region_2->get_num_cells()+1);
                        quad_children_2[i]->psi[k](0, spat_pos_2) = quad_children_1[i]->psi[k](0, spat_pos_1);
                    }
                }
            }
            // If there are no children in the temporary children vector the bottom is reached
            if (quad_children_tmp_all_1.size() == 0 && quad_children_tmp_all_2.size() == 0) {
                reached_bottom = true;
            }
            // Set the children vector as the temporary children vector for the next level calculation
            quad_children_1 = quad_children_tmp_all_1;
            quad_children_2 = quad_children_tmp_all_2;
            quad_children_tmp_all_1.clear();
            quad_children_tmp_all_2.clear();
            quad_children_tmp_ind_1.clear();
            quad_children_tmp_ind_2.clear();
        }
    }
    
    /* Case 2: Mapping from a refined quadrature to a base quadrature */
    else if (quad_selected_1->get_has_children() == true && quad_selected_2->get_has_children() == false) {
        // Set the child node as the inner root node
        quad_inner_children = quad_children_1; // Inner root node's first level children
        // Add the basis and quadrature of the root node
        basis_tmp.push_back(vector<mat>());
        wgt_tmp.push_back(vector<vec>());
        dir_tmp.push_back(vector<vec>());
        basis_tmp[0].push_back(quad_selected_1->get_basis());
        wgt_tmp[0].push_back(quad_selected_1->get_wgts());
        dir_tmp[0].push_back(quad_selected_1->get_dirs());
        // Add space to store information for the next level of children
        basis_tmp.push_back(vector<mat>());
        wgt_tmp.push_back(vector<vec>());
        dir_tmp.push_back(vector<vec>());
        // Traverse to the bottom of the inner root node
        level = 1;
        reached_bottom_inner = false;
        while (reached_bottom_inner == false) {
            // Go through the inner children nodes in the current level
            for (int k=0; k<quad_inner_children.size(); k++) {
                // If the current child has children add them to a temporary children vector
                if (quad_inner_children[k]->get_has_children() == true) {
                    quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                    // Store the basis and quadrature for mapping
                    basis_tmp[level].push_back(quad_inner_children[k]->get_basis());
                    wgt_tmp[level].push_back(quad_inner_children[k]->get_wgts());
                    dir_tmp[level].push_back(quad_inner_children[k]->get_dirs());
                }
                // Store the psi, level, column position, basis and quadrature of nodes without children
                else {
                    // Store level and column position
                    map_tmp(0,0) = level;
                    map_tmp(0,1) = k;
                    // Store psi
                    for (int r=0; r<2; r++) {
                        map_tmp(1,r) = quad_inner_children[k]->psi[r](0,spat_pos_1);
                    }
                    // Store directions
                    dirs_tmp_ind = quad_inner_children[k]->get_dirs();
                    map_tmp(2,0) = dirs_tmp_ind(0);
                    map_tmp(2,1) = dirs_tmp_ind(1);
                    // Store weights
                    dir_wgts_ind = quad_inner_children[k]->get_wgts();
                    map_tmp(3,0) = dir_wgts_ind(0);
                    map_tmp(3,1) = dir_wgts_ind(1);
                    // Store basis
                    basis_ind = quad_inner_children[k]->get_basis();
                    map_tmp(4,0) = basis_ind(0,0);
                    map_tmp(4,1) = basis_ind(0,1);
                    map_tmp(5,0) = basis_ind(1,0);
                    map_tmp(5,1) = basis_ind(1,1);
                    map_values.push_back(map_tmp);
                }
            }
            // Check if the bottom of the inner tree has been reached
            if (quad_inner_children_tmp_all.size() == 0) {
                reached_bottom_inner = true;
            }
            else {
                level++;
                basis_tmp.push_back(vector<mat>());
                wgt_tmp.push_back(vector<vec>());
                dir_tmp.push_back(vector<vec>());
            }
            // Set the inner children vector as the temporary children vetor for the next level calculation
            quad_inner_children = quad_inner_children_tmp_all;
            quad_inner_children_tmp_all.clear();
            quad_inner_children_tmp_ind.clear();
        }
        // Mapped stored psi back to child node
        // Go through all the levels from the bottom up
        for (int k=level; k>0; k--) {
            // Track coarse nodes in current level
            counter = 0;
            // Go through the map_values vector
            for (int r=0; r<map_values.size(); r++) {
                // Find the first element in map_values at the current level
                if (int(map_values[r](0,0))==k) {
                    // Create the b-vector
                    b_vec[0] = map_values[r](3,0)*map_values[r](1,0) + \
                    map_values[r](3,1)*map_values[r](1,1) + \
                    map_values[r+1](3,0)*map_values[r+1](1,0) + \
                    map_values[r+1](3,1)*map_values[r+1](1,1);
                    b_vec[1] = map_values[r](2,0)*map_values[r](3,0)*map_values[r](1,0) + \
                    map_values[r](2,1)*map_values[r](3,1)*map_values[r](1,1) + \
                    map_values[r+1](2,0)*map_values[r+1](3,0)*map_values[r+1](1,0) + \
                    map_values[r+1](2,1)*map_values[r+1](3,1)*map_values[r+1](1,1);
                    // Create the a-matrix
                    a_mat(0,0) = (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](0))*wgt_tmp[k-1][counter](0) + \
                    (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](1))*wgt_tmp[k-1][counter](1);
                    a_mat(0,1) = (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](0))*wgt_tmp[k-1][counter](0) + \
                    (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](1))*wgt_tmp[k-1][counter](1);
                    a_mat(1,0) = (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](0))*dir_tmp[k-1][counter](0)*wgt_tmp[k-1][counter](0) + \
                    (basis_tmp[k-1][counter](0,0)+basis_tmp[k-1][counter](1,0)*dir_tmp[k-1][counter](1))*dir_tmp[k-1][counter](1)*wgt_tmp[k-1][counter](1);
                    a_mat(1,1) = (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](0))*dir_tmp[k-1][counter](0)*wgt_tmp[k-1][counter](0) + \
                    (basis_tmp[k-1][counter](0,1)+basis_tmp[k-1][counter](1,1)*dir_tmp[k-1][counter](1))*dir_tmp[k-1][counter](1)*wgt_tmp[k-1][counter](1);
                    coarse_sol = solve(a_mat, b_vec);
                    // Delete the two elements that were averaged together
                    map_values.erase(map_values.begin()+r, map_values.begin()+r+2);
                    // Find col_actual which is the position in map_values to insert the new map_tmp
                    // Find col_loc which is the position relative to other elements in the same level
                    found = false;
                    at_beginning = true;
                    // Start searching for the leftmost element in the same level
                    col_loc = 0;
                    // Go through the entire map_values vector
                    for (int t=0; t<map_values.size(); t++) {
                        // Find the elements at the level of the new map_tmp
                        if (map_values[t](0,0)==k-1) {
                            // If such an element is found then col_actual cannot be at the beginning of map_values
                            at_beginning = false;
                            // Check if the element is at the current col_loc
                            if (map_values[t](0,1)==col_loc) {
                                col_loc++;        // Incrememnt col_loc
                                col_actual = t+1; // Set col_actual as the adjacent position
                                found = true;
                            }
                            // Check if the element should be inserted between two elements at the level of the new map_tmp
                            else if (map_values[t](0,1)>col_loc) {
                                col_actual = t;             // Set col_actual as the current position
                                t = int(map_values.size()); // Exit loop
                                found = true;
                            }
                        }
                    }
                    // The new map_tmp should be at the beginning of the map_values vector
                    if (at_beginning == true) {
                        col_actual = 0;
                    }
                    // The new map_tmp should be at the end of the map_values vector
                    else if (found == false) {
                        col_actual = int(map_values.size()) ;
                    }
                    // Create the new map_tmp
                    map_tmp(0,0) = k-1;
                    map_tmp(0,1) = col_loc;
                    map_tmp(1,0) = coarse_sol(0);
                    map_tmp(1,1) = coarse_sol(1);
                    map_tmp(2,0) = dir_tmp[k-1][counter](0);
                    map_tmp(2,1) = dir_tmp[k-1][counter](1);
                    map_tmp(3,0) = wgt_tmp[k-1][counter](0);
                    map_tmp(3,1) = wgt_tmp[k-1][counter](1);
                    map_tmp(4,0) = basis_tmp[k-1][counter](0,0);
                    map_tmp(4,1) = basis_tmp[k-1][counter](0,1);
                    map_tmp(5,0) = basis_tmp[k-1][counter](1,0);
                    map_tmp(5,1) = basis_tmp[k-1][counter](1,1);
                    // Insert the new map_tmp into map_values
                    map_values.insert(map_values.begin()+col_actual, map_tmp);
                    // Reset r
                    r=0;
                    // Increment counter to go to next coarse node in current level
                    counter = counter + 1;
                }
            }
        }
        map_values.clear();
        basis_tmp.clear();
        wgt_tmp.clear();
        dir_tmp.clear();
        // Pass the mapped values to the downstream LDFE region
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0,spat_pos_2) = map_tmp(1,k);
        }
        quad_inner_children.clear();
    }
    
    /* Case 3: Mapping from a base quadrature to a refined quadrature */
    else if (quad_selected_1->get_has_children() == false && quad_selected_2->get_has_children() == true) {
        // Pass solution to quadrature you are mapping to
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0, spat_pos_2) = quad_selected_1->psi[k](0, spat_pos_1);
        }
        // Set the child node as the root node
        quad_inner_children = quad_selected_2->get_children(); // Root node's first level of children
        // Map to the first level children
        basis_tmp_coarse = quad_selected_2->get_basis();
        wgts_tmp_coarse = quad_selected_2->get_wgts();
        dirs_tmp_coarse = quad_selected_2->get_dirs();
        sol_tmp_coarse = quad_selected_2->psi;
        dirs_tmp_fine_0 = quad_inner_children[0]->get_dirs();
        dirs_tmp_fine_1 = quad_inner_children[1]->get_dirs();
        wgts_tmp_fine_0 = quad_inner_children[0]->get_wgts();
        wgts_tmp_fine_1 = quad_inner_children[1]->get_wgts();
        basis_tmp_fine_0 = quad_inner_children[0]->get_basis();
        basis_tmp_fine_1 = quad_inner_children[1]->get_basis();
        // Assemble A-matrix
        a_mat(0,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
        wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
        wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
        wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1));
        a_mat(0,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))+\
        wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))+\
        wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))+\
        wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
        a_mat(1,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
        wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
        wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
        wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
        a_mat(1,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
        wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
        wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
        wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
        // Assemble b-vector
        b_vec(0) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)+\
        wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2);
        b_vec(1) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)*dirs_tmp_coarse(0)+\
        wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2)*dirs_tmp_coarse(1);
        // Solve for coarse tilde
        coarse_tilde = solve(a_mat, b_vec);
        // Solve for fine psi
        fine_sol_0(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
        coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0));
        fine_sol_0(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
        coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1));
        fine_sol_1(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
        coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0));
        fine_sol_1(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))+\
        coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
        // Store fine psi into next level children
        quad_inner_children[0]->psi.resize(2);
        quad_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[0]->psi[0](0, spat_pos_2) = fine_sol_0(0);
        quad_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[0]->psi[1](0, spat_pos_2) = fine_sol_0(1);
        quad_inner_children[1]->psi.resize(2);
        quad_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[1]->psi[0](0, spat_pos_2) = fine_sol_1(0);
        quad_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
        quad_inner_children[1]->psi[1](0, spat_pos_2) = fine_sol_1(1);
        // Traverse to the bottom of the root node
        reached_bottom_inner = false;
        while (reached_bottom_inner == false) {
            for (int k=0; k<quad_inner_children.size(); k++) {
                // If the current child has children add them to a temporary children vector
                if (quad_inner_children[k]->get_has_children() == true) {
                    quad_inner_children_tmp_ind = quad_inner_children[k]->get_children();
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[0]);
                    quad_inner_children_tmp_all.push_back(quad_inner_children_tmp_ind[1]);
                    // Add the parent flux to children nodes
                    vector<LDFE_reg*> quad_inner_inner_children = quad_inner_children[k]->get_children();
                    // Map to the first level children
                    basis_tmp_coarse = quad_inner_children[k]->get_basis();
                    wgts_tmp_coarse = quad_inner_children[k]->get_wgts();
                    dirs_tmp_coarse = quad_inner_children[k]->get_dirs();
                    sol_tmp_coarse = quad_inner_children[k]->psi;
                    dirs_tmp_fine_0 = quad_inner_inner_children[0]->get_dirs();
                    dirs_tmp_fine_1 = quad_inner_inner_children[1]->get_dirs();
                    wgts_tmp_fine_0 = quad_inner_inner_children[0]->get_wgts();
                    wgts_tmp_fine_1 = quad_inner_inner_children[1]->get_wgts();
                    basis_tmp_fine_0 = quad_inner_inner_children[0]->get_basis();
                    basis_tmp_fine_1 = quad_inner_inner_children[1]->get_basis();
                    // Assemble A-matrix
                    a_mat(0,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1));
                    a_mat(0,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                    a_mat(1,0) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                    a_mat(1,1) = wgts_tmp_fine_0(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0))*dirs_tmp_fine_0(0)+\
                    wgts_tmp_fine_0(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1))*dirs_tmp_fine_0(1)+\
                    wgts_tmp_fine_1(0)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0))*dirs_tmp_fine_1(0)+\
                    wgts_tmp_fine_1(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1))*dirs_tmp_fine_1(1);
                    // Assemble b-vector
                    b_vec(0) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)+\
                    wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2);
                    b_vec(1) = wgts_tmp_coarse(0)*sol_tmp_coarse[0](0,spat_pos_2)*dirs_tmp_coarse(0)+\
                    wgts_tmp_coarse(1)*sol_tmp_coarse[1](0,spat_pos_2)*dirs_tmp_coarse(1);
                    // Solve for coarse tilde
                    coarse_tilde = solve(a_mat, b_vec);
                    // Solve for fine psi
                    fine_sol_0(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(0))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(0));
                    fine_sol_0(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_0(1))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_0(1));
                    fine_sol_1(0) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(0))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(0));
                    fine_sol_1(1) = coarse_tilde(0)*(basis_tmp_coarse(0,0)+basis_tmp_coarse(1,0)*dirs_tmp_fine_1(1))+\
                    coarse_tilde(1)*(basis_tmp_coarse(0,1)+basis_tmp_coarse(1,1)*dirs_tmp_fine_1(1));
                    // Store fine psi into next level children
                    quad_inner_inner_children[0]->psi.resize(2);
                    quad_inner_inner_children[0]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[0]->psi[0](0, spat_pos_2) = fine_sol_0(0);
                    quad_inner_inner_children[0]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[0]->psi[1](0, spat_pos_2) = fine_sol_0(1);
                    quad_inner_inner_children[1]->psi.resize(2);
                    quad_inner_inner_children[1]->psi[0].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[1]->psi[0](0, spat_pos_2) = fine_sol_1(0);
                    quad_inner_inner_children[1]->psi[1].resize(2, region_2->get_num_cells()+1);
                    quad_inner_inner_children[1]->psi[1](0, spat_pos_2) = fine_sol_1(1);
                }
            }
            // Check if the bottom of the inner tree has been reached
            if (quad_inner_children_tmp_all.size() == 0) {
                reached_bottom_inner = true;
            }
            // Set the inner children vector as the temporary children vetor for the next calculation
            quad_inner_children = quad_inner_children_tmp_all;
            quad_inner_children_tmp_all.clear();
            quad_inner_children_tmp_ind.clear();
        }
        quad_inner_children.clear();
    }
    
    /* Case 4: Mapping between two base quadratures */
    else {
        // Pass the solution to the quadrature you are mapping to
        quad_selected_2->psi.resize(2);
        for (int k=0; k<2; k++) {
            quad_selected_2->psi[k].resize(2, region_2->get_num_cells()+1);
            quad_selected_2->psi[k](0, spat_pos_2) = quad_selected_1->psi[k](0, spat_pos_1);
        }
    }
}
