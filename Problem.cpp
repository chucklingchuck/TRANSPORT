//
//  Problem.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/15/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Problem.h"

Problem::Problem (Input input_param) {
    
    /* Retrieve input parameters */
    int num_reg         = input_param.get_num_reg();
    string quad_type    = input_param.get_quad_type();
    Col<int> num_cells  = input_param.get_num_cells();
    Col<int> quad_order = input_param.get_quad_order();
    Col<int> scat_order = input_param.get_scat_order();
    vec reg_size        = input_param.get_reg_size();
    vec abs_xs          = input_param.get_abs_xs();
    vec scat_xs_all     = input_param.get_scat_xs();
    vec ext_source      = input_param.get_ext_source();
    vector<vec> scat_xs_ind;
    int counter = 0;
    for (int i=0; i<num_reg; i++){
        scat_xs_ind.push_back(vec());
        scat_xs_ind[i].resize(scat_order[i]+1);
        for (int j=0; j<scat_order[i]+1; j++) {
            scat_xs_ind[i](j)=scat_xs_all(counter);
            counter++;
        }
    }
    
    /* Create vector of pointers to each region */
    Regions.resize(num_reg);
    for (int i=0; i<num_reg; i++){
        Regions[i]=new Region(quad_order(i), num_cells(i), reg_size(i), scat_order(i), \
                              abs_xs(i), scat_xs_ind[i], ext_source(i), quad_type);
    }

}
