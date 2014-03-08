//
//  Output.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/25/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Output.h"

void Output(Problem problem) {

    /* Define basic parameters used */
    vector<Region*> Regions = problem.get_Regions();
    int num_reg = int(Regions.size());
    int num_cells;
    int num_edges;
    int counter;
    char buffer[32];
    bool reached_bottom;
    LDFE_quad* quad_region;
    LDFE_reg* quad_selected;
    cube psi_store;
    vec reg_size_store(num_reg);
    vec abs_xs_store(num_reg);
    vec ext_source_store(num_reg);
    mat wgts_store;
    mat dirs_store;
    Col<int> num_cells_store(num_reg);
    vector<LDFE_reg*> quad_ind;
    vector<LDFE_reg*> quad_children;
    vector<LDFE_reg*> quad_children_tmp_ind;
    vector<LDFE_reg*> quad_children_tmp_all;
    vector<mat> psi_store_tmp;
    vector<vec> wgts_store_tmp;
    vector<vec> dirs_store_tmp;
    
    /* Write psi, weights and directions of each region to temporary vectors */
    for (int i=0; i<num_reg; i++) {
        num_cells = Regions[i]->get_num_cells();
        num_edges = num_cells+1;
        for (int j=0; j<2; j++) {
            quad_region   = Regions[i]->get_quad();    // Region quadrature pointer
            quad_ind      = quad_region->get_quad();   // Vector of pos and neg LDFE regions
            quad_selected = quad_ind[j];               // Selected root LDFE node
            if (quad_selected->get_has_children()==true) {
                quad_children = quad_selected->get_children(); // Root node's 1st level children
                // Traverse the tree
                reached_bottom = false;
                while (reached_bottom == false) {
                    // Go through all the children nodes of current level
                    for (int j=0; j<quad_children.size(); j++) {
                        // If the child has children add them to a temporary children vector
                        if (quad_children[j]->get_has_children() == true) {
                            quad_children_tmp_ind = quad_children[j]->get_children();
                            quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
                            quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
                        }
                        // If child does not have children store the data in the temporary holders
                        else {
                            psi_store_tmp.push_back(quad_children[j]->psi[0]);
                            psi_store_tmp.push_back(quad_children[j]->psi[1]);
                            wgts_store_tmp.push_back(quad_children[j]->get_wgts());
                            dirs_store_tmp.push_back(quad_children[j]->get_dirs());
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
            // Directly store the data into temporary holder for the coarsest quadrature
            else {
                psi_store_tmp.push_back(quad_selected->psi[0]);
                psi_store_tmp.push_back(quad_selected->psi[1]);
                wgts_store_tmp.push_back(quad_selected->get_wgts());
                dirs_store_tmp.push_back(quad_selected->get_dirs());
            }
        }
        
        /* Convert to arma format for HDF5 storage */
        psi_store.resize(2, num_edges, int(psi_store_tmp.size()));
        wgts_store.resize(int(wgts_store_tmp.size()), 2);
        dirs_store.resize(int(dirs_store_tmp.size()), 2);
        counter = 0;
        for (int j=0; j<wgts_store_tmp.size(); j++){
            psi_store.slice(counter) = psi_store_tmp[counter];
            psi_store.slice(counter+1) = psi_store_tmp[counter+1];
            for (int k=0; k<2; k++) {
                wgts_store(j, k) = wgts_store_tmp[j](k);
                dirs_store(j, k) = dirs_store_tmp[j](k);
            }
            counter = counter + 2;
        }
        psi_store_tmp.clear();
        dirs_store_tmp.clear();
        wgts_store_tmp.clear();
        snprintf(buffer, sizeof(char) * 32, "psi_reg%i", i);
        psi_store.save(buffer, hdf5_binary);
        snprintf(buffer, sizeof(char) * 32, "dirs_reg%i", i);
        dirs_store.save(buffer, hdf5_binary);
        snprintf(buffer, sizeof(char) * 32, "wgts_reg%i", i);
        wgts_store.save(buffer, hdf5_binary);
        psi_store.clear();
        dirs_store.clear();
        wgts_store.clear();
    }
    
    /* Write other problem parameters to HDF5 files */
    for (int i=0; i<num_reg; i++) {
        reg_size_store(i)   = Regions[i]->get_reg_size();
        num_cells_store(i)  = Regions[i]->get_num_cells();
        abs_xs_store(i)     = Regions[i]->get_abs_xs();
        ext_source_store(i) = Regions[i]->get_ext_source();
    }
    reg_size_store.save("reg_size", hdf5_binary);
    num_cells_store.save("num_cells", hdf5_binary);
    abs_xs_store.save("abs_xs",hdf5_binary);
    ext_source_store.save("ext_source",hdf5_binary);
    
}
