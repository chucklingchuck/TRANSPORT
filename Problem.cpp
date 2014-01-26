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
    int num_reg = input_param.get_num_reg();
    Col<int> num_cells = input_param.get_num_cells();
    vec reg_size = input_param.get_reg_size();
    string quad_type = input_param.get_quad_type();
    Col<int> quad_order = input_param.get_quad_order();
    vec abs_xs = input_param.get_abs_xs();
    vec scat_xs = input_param.get_scat_xs();
    vec ext_source = input_param.get_ext_source();
    
    /* Create vector of pointers to each region */
    Regions.resize(num_reg);
    for (int i=0; i<num_reg; i++){
        Regions[i] = new Region(quad_order(i), num_cells(i), reg_size(i), \
                                abs_xs(i), scat_xs(i), ext_source(i), quad_type);
    }
    
}
