//
//  Quadrature.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/11/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//
//  INFORMATION: DOCUMENTx ATION/ANGULAR/

#include "Quadrature.h"

LDFE_quad::LDFE_quad(int param) {
    
    /* Calculate number of LDFE ranges */
    order = param;
    int num_LDFE = 2 * pow(2, order);
    
    /* Calculate quadrature directions */
    double delta_mu = (2./num_LDFE)/4.;
    double mu_current = -1.;
    mat mu(num_LDFE, 2);
    for (int i=0; i<num_LDFE; i++) {
        mu(i, 0) = mu_current+delta_mu;
        mu(i, 1) = mu(i, 0)+2.*delta_mu;
        mu_current = mu_current+4.*delta_mu;
    }
    
    /* Solve the basis functions of each LDFE range */
    mat A = mat(2, 2);
    mat B = eye<mat>(2, 2);
    mat C = mat(3, 3);
    vec wgt = vec(2);
    mat dirs_temp = mat(num_LDFE, 2);
    mat wgts_temp = mat(num_LDFE, 2);
    wgts.set_size(num_LDFE, 2);
    dirs.set_size(num_LDFE, 2);
    basis.resize(num_LDFE);
    for (int i=0; i<num_LDFE; i++){
        for (int j=0; j<2; j++){
            A(j, 0) = 1.;
            A(j, 1) = mu(i, j);
        }
        C = solve(A, B);
        
        /* Solve for the weights of each LDFE range */
        double mu_max = mu(i, 1)+delta_mu;
        double mu_min = mu(i, 0)-delta_mu;
        for (int j=0; j<2; j++){
            wgt(j) = C(0, j)*(mu_max-mu_min)+ \
            C(1, j)*(pow(mu_max,2.)-pow(mu_min, 2))/2.;
        }
        
        /* Store quadrature directions, weights and basis functions */
        for (int j=0; j<2; j++) {
            dirs_temp(i, j) = mu(i, j);
            wgts_temp(i, j) = wgt(j);
        }
        basis[i] = C;
    }
    
    /* Flip weights and direction matrices */
    int counter = num_LDFE-1;
    for (int i=0; i<num_LDFE; i++) {
        dirs(i,0) = dirs_temp(counter,1);
        dirs(i,1) = dirs_temp(counter,0);
        wgts(i,0) = wgts_temp(counter,1);
        wgts(i,1) = wgts_temp(counter,0);
        counter = counter - 1;
    }
}

