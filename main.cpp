//
//  main.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/11/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Common.h"
#include "Input.h"
#include "Problem.h"
#include "Solver.h"
#include "Output.h"

int main()
{
    
    /* Read in XML input file */
    cout << "Reading in XML input file..." << endl;
    Input input_param;

    /* Build problem */
    cout << "Building the problem..." << endl;
    Problem problem(input_param);
        
    /* Call solver */
    cout << "Running the solver..." << endl;
    Solver(input_param, problem);

    /* Write HDF5 output file */
    cout << "Writing the output..." << endl;
    Output(problem);

    /* Exit program */
    cout << "Program successfully ran!" << endl;
    EXIT_SUCCESS;
    
}

