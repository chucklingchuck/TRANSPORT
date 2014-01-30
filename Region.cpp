//
//  Region.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Region.h"

Region::Region (int param1, int param2, double param3, int param4, double param5, vec param6, double param7, string param8) {
    
    /* Store region parameters */
    order = param1;
    num_cells = param2;
    reg_size = param3;
    scat_order = param4;
    abs_xs = param5;
    scat_xs = param6;
    ext_source = param7;
    quad_type = param8;
    
    /* Create quadrature */
    if (quad_type == "LDFE") {
        quad_ptr = new LDFE_quad(order);
    }
    else {
        exit(EXIT_FAILURE);
    }
    
}
