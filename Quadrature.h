//
//  Quadrature.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/11/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef __TRANSPORT__Quadrature__
#define __TRANSPORT__Quadrature__

#include "Common.h"

class Quadrature {
protected:
    mat dirs;
    mat wgts;
    int order;
    vector<mat> basis;
public:
    mat get_dirs() {return dirs;};
    mat get_wgts() {return wgts;};
    vector<mat> get_basis() {return basis;};
};

class LDFE_quad: public Quadrature {
public:
    LDFE_quad(int);
};

#endif /* defined(__TRANSPORT__Quadrature__) */
