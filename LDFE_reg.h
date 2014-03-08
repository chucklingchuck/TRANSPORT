//
//  LDFE_reg.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 2/1/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_LDFE_reg_h
#define TRANSPORT_LDFE_reg_h

#include "Common.h"

class LDFE_reg {
protected:
    vec dirs;
    vec wgts;
    mat basis;
public:
    bool has_children;
    vector<LDFE_reg*> children;
    bool get_has_children() {return has_children;};
    vector<LDFE_reg*> get_children() {return children;};
    vec get_dirs() {return dirs;};
    vec get_wgts() {return wgts;};
    mat get_basis() {return basis;};
    vector<mat> psi;
    LDFE_reg(double, double, int, int, string);
};

#endif
