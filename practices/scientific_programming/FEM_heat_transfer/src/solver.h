#ifndef SOLVER_H_
#define SOLVER_H_

#include "node.h"
#include "element.h"
#include "postProcessor.h"
#include "dco.hpp"
using namespace dco;
typedef gt1s<double>::type DCO_TYPE;


/***************************************************************************************************
 * Solver class
***************************************************************************************************/
class femSolver
{
    private:
        setting*    stngs;  // a local pointer for the settings
        node*       nodes;  // a local pointer for the nodes
        element*    elems;  // a local pointer for the elements
        postProcessor postP;   // Initialize postprocessor
    // AD variables
      DCO_TYPE* Ta;
      DCO_TYPE ka;
      DCO_TYPE va;

        // PRIVATE METHODS
        void explicitSolver();

    protected:

    public:
        // CONSTRUCTOR //
        femSolver(setting* argSet, node* argNode, element* argElem)
        {
            stngs = argSet;
            nodes = argNode;
            elems = argElem;
            Ta = new DCO_TYPE[stngs->getNn()];
            
        };

        // DESTRUCTOR
        ~femSolver(){
        delete[] Ta;
        };
        
        // INTERFACE FUNCTION
        void solve();

};

#endif /* SOLVER_H_ */


