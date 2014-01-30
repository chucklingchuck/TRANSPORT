//
//  Sweep.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Sweep_h
#define TRANSPORT_Sweep_h

#include "Common.h"
#include "Region.h"

class Sweep {
protected:
    vector<mat> region_psi;
public:
    vector<mat> get_region_psi() {return region_psi;};
};

class LDFE_spat: public Sweep {
public:
    LDFE_spat(Region*, vector<mat>, vector<mat>, string);
};

double legendre (double, double);

#endif
