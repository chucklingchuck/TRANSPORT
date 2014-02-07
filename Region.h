//
//  Region.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Region_h
#define TRANSPORT_Region_h

#include "Common.h"
#include "Quadrature.h"

class Region {
protected:
    int order;
    int num_cells;
    int scat_order;
    double reg_size;
    double abs_xs;
    double ext_source;
    string quad_type;
    LDFE_quad* quad_ptr;
    vec scat_xs;
public:
    int get_order() {return order;};
    int get_num_cells() {return num_cells;};
    int get_scat_order() {return scat_order;};
    double get_reg_size() {return reg_size;};
    double get_abs_xs() {return abs_xs;};
    double get_ext_source() {return ext_source;};
    string get_quad_type() {return quad_type;};
    LDFE_quad* get_quad() {return quad_ptr;};
    vec get_scat_xs() {return scat_xs;};
    Region (int, int, double, int, double, vec, double, string);
};

#endif
