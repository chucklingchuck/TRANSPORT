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

class Mapping {
protected:
    vec psi_to;
public:
    vec get_psi_to() {return psi_to;};
};

class one_to_one: public Mapping {
public:
    one_to_one(vector<mat>, string);
};

/*
class Jarrel1: public Mapping {
public:
    Jarrel1(vector<mat>, int, int);
};
*/

#endif
