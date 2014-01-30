//
//  main.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/11/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Input.h"
#include "Output.h"
#include "Problem.h"
#include "Solver.h"
#include "Common.h"

int main()
{
    
    /* Read in XML input file */
    Input input_param;

    /* Build problem */
    Problem problem_setup(input_param);

    /* Call solver */
    Solver solution(input_param, problem_setup);

    /* Write HDF5 output file */
    Output(problem_setup, solution);

    /* Exit program */
    EXIT_SUCCESS;
    
}

