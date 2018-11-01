/***************************************************************************************************
 * Name        : 2D_Unsteady_Diffusion.cpp
 * Author      : A. Emre Ongut
 * Version     : 2.0
 * Copyright   : See the copyright notice in the README file.
 * Description : This code is designed for educational purposes.
 *               - For now, everything is for 2D linear triangular elements
 *               - uses Gauss quadrature integration with 7 points
 *               - supports constant Drichlet and Neumann boundary conditions
 *               - not yet in parallel
***************************************************************************************************/

#include "setting.h"
#include "node.h"
#include "element.h"
#include "solver.h"
#include "postProcessor.h"
// I added more
#include "cmath"
#include "cstdlib"
#include "string"
#include "vector"
#include <iostream>

int main(void)
{
    /***********************************************************************************************
     *  Pre-processing stage: 
     * Read input files, allocate and initialize arrays, etc.
     **********************************************************************************************/
    setting stngs;              // A setting object is createad.
    stngs.readSettings();       // settings.in file and minf file are read.

    node nodes(stngs.getNn()); // We create a node object
    nodes.createNodes(&stngs); // and read nodal data from file.

    element elems(stngs.getNe(), nenTri, nefTri, 7);     // We create an element object. (,,7 point gaussq)
    elems.createElements(&stngs, &nodes);                // and read element data from file
    // Moreover solution independent calculations are made: Jacobian determinant, shape funct.
    // derivatives as well as element matrices are calculated.


    /***********************************************************************************************
     * Solution stage: 
     * Element matrices are calculated, boundary conditions are applied and system is solved.
     **********************************************************************************************/

    femSolver solver(&stngs, &nodes, &elems);   // We create an femSolver object.
    solver.solve();                             // and solve the equation system

    /***********************************************************************************************
     * Post-Processing Stage
     * Create the output files / visualize the reuslts.
     **********************************************************************************************/




    cout << endl << "Ciao :)" << endl;
    return 0;
}

