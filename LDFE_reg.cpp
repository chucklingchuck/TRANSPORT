//
//  LDFE_reg.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 2/1/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "LDFE_reg.h"

LDFE_reg::LDFE_reg (double min_mu, double max_mu, int current_order, int final_order) {

    /* Calculate quadrature directions */
    dirs.resize(2);
    wgts.resize(2);
    children.resize(2);
    basis.resize(2,2);
    dirs(0) = min_mu+(max_mu-min_mu)/4.;
    dirs(1) = max_mu-(max_mu-min_mu)/4.;
    
    /* Solve the basis functions of each LDFE range */
    mat A = mat(2, 2);
    mat B = eye<mat>(2, 2);
    basis.resize(3, 3);
    for (int j=0; j<2; j++){
        A(j, 0) = 1.;
        A(j, 1) = dirs(j);
    }
    basis = solve(A, B);
    
    /* Solve for the weights of each LDFE range */
    for (int j=0; j<2; j++){
        wgts(j) = basis(0, j)*(max_mu-min_mu)+ \
                  basis(1, j)*(pow(max_mu,2.)-pow(min_mu, 2))/2.;
    }
    
    /* Generate children LDFE regions */
    if (current_order == final_order) {
        has_children = false;
        children[0] = NULL;
        children[1] = NULL;
    }
    else {
        has_children = true;
        double mu_center = (max_mu+min_mu)/2.;
        children[0] = new LDFE_reg(min_mu, mu_center, current_order+1, final_order); // Left child
        children[1] = new LDFE_reg(mu_center, max_mu, current_order+1, final_order); // Right child
    }
}