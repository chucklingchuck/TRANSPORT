//
//  Region.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Region_h
#define TRANSPORT_Region_h

#include "Quadrature.h"
#include "Common.h"

class Region {
protected:
    Quadrature* quad_ptr;
    int order;
    int num_cells;
    double reg_size;
    double abs_xs;
    double scat_xs;
    double ext_source;
    string quad_type;
public:
    int get_order() {return order;};
    int get_num_cells() {return num_cells;};
    double get_reg_size() {return reg_size;};
    double get_abs_xs() {return abs_xs;};
    double get_scat_xs() {return scat_xs;};
    double get_ext_source() {return ext_source;};
    string get_quad_type() {return quad_type;};
    mat get_quad_wgts() {return quad_ptr->get_wgts();};
    mat get_quad_dirs() {return quad_ptr->get_dirs();};
    vector<mat> get_quad_basis() {return quad_ptr->get_basis();};
    Region(int, int, double, double, double, double, string);
};

#endif
