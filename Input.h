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
    int num_reg;
    int si_cycles;
    int max_cycles_per_adapt;
    int si_per_adapt;
    double si_tol;
    double psi_left;
    double psi_right;
    double alpha;
    double beta;
    double adapt_tol;
    double eps_thresh;
    string spat_type;
    string map_type;
    string fixup_type;
    string quad_type;
    string bc_type;
    string adapt_type;
    Col<int> num_cells;
    Col<int> quad_order;
    Col<int> scat_order;
    vec abs_xs;
    vec scat_xs;
    vec reg_size;
    vec ext_source;
public:
    int get_num_reg() {return num_reg;};
    int get_si_cycles() {return si_cycles;};
    int get_si_per_adapt() {return si_per_adapt;};
    int get_max_cycles_per_adapt() {return max_cycles_per_adapt;};
    double get_si_tol() {return si_tol;};
    double get_psi_left() {return psi_left;};
    double get_psi_right() {return psi_right;};
    double get_alpha() {return alpha;};
    double get_beta() {return beta;};
    double get_adapt_tol() {return adapt_tol;};
    double get_eps_thresh() {return eps_thresh;};
    string get_spat_type() {return spat_type;};
    string get_map_type() {return map_type;};
    string get_fixup_type() {return fixup_type;};
    string get_quad_type() {return quad_type;};
    string get_bc_type() {return bc_type;};
    string get_adapt_type() {return adapt_type;};
    Col<int> get_num_cells() {return num_cells;};
    Col<int> get_quad_order() {return quad_order;};
    Col<int> get_scat_order() {return scat_order;};
    vec get_abs_xs() {return abs_xs;};
    vec get_scat_xs() {return scat_xs;};
    vec get_reg_size() {return reg_size;};
    vec get_ext_source() {return ext_source;};
    Input();
};

#endif
