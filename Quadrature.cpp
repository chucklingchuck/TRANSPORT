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
    
    /* Calculate number of LDFE ranges */
    order = param;                    // Quadrature order
    int num_LDFE = 2 * pow(2, order); // Number of LDFE regions
    
    /* Calculate quadrature directions */
    double delta_gamma = (pi / num_LDFE) / 4.0;            // Gamma from node to nearest LDFE edge
    double gamma_current = 0.0;
    mat gamma(num_LDFE, 2);
    for (int i=0; i<num_LDFE; i++) {
        gamma(i, 0) = gamma_current + delta_gamma;         // Left node gamma
        gamma(i, 1) = gamma(i, 0) + 2.0 * delta_gamma;     // Right node gamma
        gamma_current = gamma_current + 4.0 * delta_gamma;
    }
    
    /* Solve the basis functions of each LDFE range */
    mat A = mat(2, 2);
    mat B = eye<mat>(2, 2);
    mat C = mat(3, 3);
    vec wgt = vec(2);
    wgts.set_size(num_LDFE, 2);          // Resize matrix storing all LDFE range weights
    dirs.set_size(num_LDFE, 2);          // Resize matrix storing all LDFE range directions
    basis.resize(num_LDFE);              // Resize vector storing all LDFE basis functions
    for (int i=0; i<num_LDFE; i++){
        for (int j=0; j<2; j++){
            A(j, 0) = 1.0;
            A(j, 1) = cos(gamma(i, j));
        }
        C = solve(A, B);                 // Solve for the basis functions

        /* Solve for the weights of each LDFE range */
        double gamma_max = gamma(i, 1) + delta_gamma;
        double gamma_min = gamma(i, 0) - delta_gamma;
        for (int j=0; j<2; j++){
            wgt(j) = C(0, j) * (gamma_max - gamma_min) + \
                     C(1, j) * (sin(gamma_max) - sin(gamma_min));
        }

        /* Store quadrature directions, weights and basis functions */
        for (int j=0; j<2; j++) {
            dirs(i, j) = cos(gamma(i, j));
            wgts(i, j) = wgt(j) / pi * 2.0;
        }
        basis[i] = C;
    }
}

LDFE_quad_eq::LDFE_quad_eq(int param) {
    
    /* Calculate number of LDFE ranges */
    order = param;                    // Quadrature order
    int num_LDFE = 2 * pow(2, order); // Number of LDFE regions
    
    /* Calculate even quadrature weights */
    double delta_gamma = pi / num_LDFE;
    double wgt_even = delta_gamma / 2.0;
    
    /* Initialize common parameters */
    double gamma_current = 0.0;
    double gamma_min;
    double gamma_max;
    double gamma_center;
    vec dir = vec(2);
    mat A = mat(2, 2);
    mat B = eye<mat>(2, 2);
    mat C = mat(3, 3);
    vec wgt = vec(2);
    wgts.set_size(num_LDFE, 2);   // Resize matrix storing all LDFE range weights
    dirs.set_size(num_LDFE, 2);   // Resize matrix storing all LDFE range directions
    basis.resize(num_LDFE);       // Resize vector storing all LDFE basis functions
    double RE_old;
    double RE_new;
    double ratio_old;
    double ratio_new;
    double ratio_temp;
    int converged;
    int iter_counter;
    double delta;
    for (int r=0; r<num_LDFE; r++) {
        
        /* First quadrature direction guess */
        gamma_min = gamma_current;
        gamma_max = gamma_current + delta_gamma;
        gamma_current = gamma_max;
        gamma_center = (gamma_max + gamma_min) / 2.0;
        ratio_old = 0.5;
        dir(0) = gamma_center - ratio_old * delta_gamma / 2.0;
        dir(1) = gamma_center + ratio_old * delta_gamma / 2.0;
        
        /* Solve for basis functions */
        for (int j=0; j<2; j++){
            A(j, 0) = 1.0;
            A(j, 1) = cos(dir(j));
        }
        C = solve(A, B); // Solve for the basis functions
        
        /* Solve for weights */
        for (int j=0; j<2; j++){
            wgt(j) = C(0, j) * (gamma_max - gamma_min) + \
                     C(1, j) * (sin(gamma_max) - sin(gamma_min));
        }
        
        /* Calculate relative error using first direction guess */
        RE_old = (wgt(0) - wgt_even) / wgt_even;
        
        /* Second ratio guess */
        ratio_new = 0.75;
        
        /* Iterate until weights are equal */
        converged = 0;
        iter_counter = 0;
        while (converged == 0 && iter_counter < 100) {
            
            /* New quadrature directions */
            dir(0) = gamma_center - ratio_new * delta_gamma / 2.0;
            dir(1) = gamma_center + ratio_new * delta_gamma / 2.0;
            
            /* New quadrature basis functions */
            for (int j=0; j<2; j++){
                A(j, 0) = 1.0;
                A(j, 1) = cos(dir(j));
            }
            C = solve(A, B); // Solve for the basis functions
            
            /* New quadrature weights */
            for (int j=0; j<2; j++){
                wgt(j) = C(0, j) * (gamma_max - gamma_min) + \
                         C(1, j) * (sin(gamma_max) - sin(gamma_min));
            }
            
            /* Calculate new relative error */
            RE_new = (wgt(0) - wgt_even) / wgt_even;
            
            /* Calculate next direction guess */
            if (RE_old == RE_new && iter_counter == 1) {
                ratio_temp = ratio_new;
                ratio_new = (ratio_old + ratio_new) / 2.0;
                ratio_old = ratio_temp;
                iter_counter++;
            }
            else {
                delta = ((ratio_new - ratio_old) / (RE_new - RE_old)) * RE_new;
                
                /* Check for convergence */
                if (abs(delta) < 1.e-12) {
                    converged = 1;
                }
                else {
                    RE_old = RE_new;
                    ratio_old = ratio_new;
                    ratio_new = ratio_new - delta;
                    iter_counter++;
                }
            }
        }
        
        /* Store quadrature directions, weights and basis functions */
        for (int j=0; j<2; j++) {
            dirs(r, j) = cos(dir(j));
            wgts(r, j) = wgt(j) / pi * 2.0;
        }
        basis[r] = C;
    }
}
