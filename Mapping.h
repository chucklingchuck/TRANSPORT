//
//  Mapping.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Mapping_h
#define TRANSPORT_Mapping_h

#include "Common.h"
#include "Region.h"

// Simple testing algorithm
void mapping1 (Region*, Region*, string);

// Preserve zeroth and first moments across LDFE regions
void mapping2 (Region*, Region*, string);

#endif
