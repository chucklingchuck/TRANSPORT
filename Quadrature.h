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
#include "LDFE_reg.h"

class LDFE_quad {
protected:
    vector<LDFE_reg*> quad;    
public:
    vector<LDFE_reg*> get_quad() {return quad;};
    LDFE_quad (int);
};

#endif
