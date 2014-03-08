//
//  Fixup.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 2/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Fixup_h
#define TRANSPORT_Fixup_h

#include "Common.h"
#include "Region.h"
#include "Sweep.h"

// Linear interpolation
void fixup1 (int, int, string, double, vector<LDFE_reg*>, vector<LDFE_reg*>);

// Left and right total source
void fixup_source (string, vec, int, vector<Region*>, vec&, vec&, vector<vector<mat> >, Col<int>, Col<int>, vec);

#endif
