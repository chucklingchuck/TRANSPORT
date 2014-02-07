//
//  Quadrature.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/11/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//
//  INFORMATION: DOCUMENTATION/ANGULAR/

#include "Quadrature.h"

LDFE_quad::LDFE_quad(int param) {
    
    /* Base quadrature parameters */
    int order = param;
    double mu_min    = -1.;
    double mu_center = 0.;
    double mu_max    = 1.;
    
    /* Recursively create LDFE quadrature */
    quad.resize(2);
    quad[0] = new LDFE_reg(mu_center, mu_max, 0, order); // Positive LDFE regions
    quad[1] = new LDFE_reg(mu_min, mu_center, 0, order); // Negative LDFE regions

}

