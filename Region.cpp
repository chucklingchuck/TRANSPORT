//
//  Region.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Region.h"

Region::Region (int param1, int param2, double param3, double param4, double param5, double param6, string param7) {
    
    /* Store region parameters */
    order = param1;
    num_cells = param2;
    reg_size = param3;
    abs_xs = param4;
    scat_xs = param5;
    ext_source = param6;
    quad_type = param7;
    
    /* Create quadrature */
    if (quad_type == "LDFE") {
        quad_ptr = new LDFE_quad(order);
    }
    else if (quad_type == "LDFE_EQ") {
        quad_ptr = new LDFE_quad_eq(order);
    }
    
}
