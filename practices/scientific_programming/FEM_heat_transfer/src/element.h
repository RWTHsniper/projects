#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "node.h"
#include "gaussq.h"

/***************************************************************************************************
 * ELEMENT LEVEL DATA STRUCTURE. 
 * Derived from GaussQ class
 ***************************************************************************************************
 * Each element has
 *      connectivity [ne x nen]
 *      face groups defined for each face [ne x nef]
 *      Mass matrix [ne x nen]
 *      Forcing matrix [ne x nen]
 *      Stiffness matrix [ne x nen x nen]
***************************************************************************************************/
class element : public GaussQ
{
    private:
        // PRIVATE VARIABLES
        size_t  ne;     // Number of elements
        size_t  nen;    // Number of element nodes
        size_t  nef;    // Number of element faces
        size_t  nGQP;   // Number of Gauss quadrature points
        int*    conn;   // Connectivity             [ne x nen]
        int*    FG;     // Face groups              [ne x nef]
        double* M;      // Element mass matrix      [ne * nen]
        double* F;      // Element RHS              [ne * nen]
        double* K;      // Element stiffness matrix [ne * nen * nen]

        // PRIVATE METHODS
        void setM (size_t ie, size_t ien, double value) {M[ie*nen+ien] = value;};	// M[ie][ien]
        void setF (size_t ie, size_t ien, double value) {F[ie*nen+ien] = value;};
        void addF (size_t ie, size_t ien, double value) {F[ie*nen+ien] += value;};
        void setK (size_t ie, size_t ien, size_t jen, double value) {K[ie*nen*nen+ien*nen+jen] = value;};

        int getFG   (size_t ie, size_t ief) {return FG[ie*nef+ief];}

        void readMien(setting*);
        void readMrng(setting*);
        void calculateJacobian(node*);
        void calculateElementMatrices(setting*, node*);
        void applyBoundaryConditions(setting*, node*);

    protected:

    public:
        // CONSTRUCTOR //	
        element(size_t argNe, size_t argNen, size_t argNef, size_t argNGQP)
        : GaussQ(argNe, argNen, argNGQP)
        {
            ne   = argNe;
            nen  = argNen;
            nef  = argNef;
            nGQP = argNGQP;
            conn = new int[ne*nen]; 
            FG   = new int[ne*nef];
            M    = new double[ne*nen];
            F    = new double[ne*nen];
            K    = new double[ne*nen*nen];
        };

        // DESTRUCTOR //
        ~element()
        {
            delete[] conn;
            delete[] FG;
            delete[] K;
            delete[] F;
            delete[] M;
        };

        // PUBLIC GETTERS
        size_t  getNen  ()          {return nen;};
        size_t  getNGQP ()          {return nGQP;};
        int     getConn (size_t ie, size_t ien) {return conn[ie*nen+ien];};
        double* getMptr (size_t ie) {return &M[ie*nen];};
        double* getFptr (size_t ie) {return &F[ie*nen];};
        double* getKptr (size_t ie) {return &K[ie*nen*nen];};

        // PUBLIC INTERFACE
        void createElements(setting*, node*);
};

#endif /* ELEMENT_H_ */
