//
//  Input.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/20/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Input_h
#define TRANSPORT_Input_h

#include "Common.h"
#include "rapidxml.hpp"

using namespace rapidxml;

class Input {
protected:
    string spat_type;
    int num_reg;
    Col<int> num_cells;
    vec reg_size;
    string quad_type;
    Col<int> quad_order;
    vec abs_xs;
    Col<int> scat_order;
    vec scat_xs;
    double si_tol;
    int si_cycles;
    vec ext_source;
    double psi_left;
    double psi_right;
public:
    int get_num_reg() {return num_reg;};
    Col<int> get_num_cells() {return num_cells;};
    vec get_reg_size() {return reg_size;};
    string get_spat_type() {return spat_type;};
    string get_quad_type() {return quad_type;};
    Col<int> get_quad_order() {return quad_order;};
    vec get_abs_xs() {return abs_xs;};
    Col<int> get_scat_order() {return scat_order;};
    vec get_scat_xs() {return scat_xs;};
    double get_si_tol() {return si_tol;};
    int get_si_cycles() {return si_cycles;};
    vec get_ext_source() {return ext_source;};
    double get_psi_left() {return psi_left;};
    double get_psi_right() {return psi_right;};
    Input();
};

#endif
