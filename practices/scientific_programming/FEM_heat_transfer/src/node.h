#ifndef NODE_H_
#define NODE_H_

#include "setting.h"

/***************************************************************************************************
 * GENERAL NODE LEVEL DATA STRUCTURE
 ***************************************************************************************************
 * Here we store the nodal data in a structure of arrays form.
 **************************************************************************************************/
class node
{
    private:
        // PRIVATE VARIABLES
        size_t  nn;     // Number of nodes
        int*    BCtype; // Boundary condition type
        double* x;      // x-coordinate array pointer
        double* y;      // y-ccordinate array pointer
        double* T;      // Temperature array pointer
        double* mass;   // nodal mass

        // PRIVATE METHODS
        void readMxyz(setting* settings);
        void readData(setting* settings);

    protected:

    public:
        // CONSTRUCTOR
        node(size_t argNn)
        {
            nn     = argNn;
            BCtype = new int[nn];
            x      = new double[nn];
            y      = new double[nn];
            T      = new double[nn];
            mass   = new double[nn];
        };

        // DESTRUCTOR
        ~node()
        {
            delete[] BCtype;
            delete[] x;
            delete[] y;
            delete[] T;
            delete[] mass;
        };

        // GETTERS
        inline double  getX(size_t i)      {return x[i];};
        inline double  getY(size_t i)      {return y[i];};
        inline double  getT(size_t i)      {return T[i];};
        inline double  getMass(size_t i)   {return mass[i];};
        inline int     getBCtype(size_t i) {return BCtype[i];};

        // SETTERS 
        inline void setMass   (size_t i, double value) {mass[i]    = value;};
        inline void addMass   (size_t i, double value) {mass[i]   += value;};
        inline void setT      (size_t i, double value) {T[i]       = value;};
        inline void setBCtype (size_t i, int    value) {BCtype[i]  = value;};

        // INTERFACE
        void createNodes(setting* settings);
};

#endif /* NODE_H_ */
