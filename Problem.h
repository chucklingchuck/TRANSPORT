//
//  Problem.h
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#ifndef TRANSPORT_Problem_h
#define TRANSPORT_Problem_h

#include "Input.h"
#include "Region.h"
#include "Common.h"

class Problem {
protected:
    vector<Region*> Regions;
public:
    vector<Region*> get_Regions() {return Regions;};
    Problem(Input);
};

#endif
