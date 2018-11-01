#ifndef SETTING_H_
#define SETTING_H_

#include "constants.h"

/***************************************************************************************************
 * Boundary Condition specification class
***************************************************************************************************/
class bndc  // size of BC is new bndc [nFG+1]
{
    private:
        friend class setting;
        int BCType;
        double BCValue1;
        double BCValue2;
    protected:

    public:
        //DEFAULT CONSTRUCTOR
        bndc(){BCType=0; BCValue1=0.0; BCValue2 = 0.0;}; // Type 1 = Drichlet, Type 2 = Neumann, Type 3 = Mixed
        //SETTERS
        void setType    (int value)     {BCType = value;};
        void setValue1  (double value)  {BCValue1 = value;};
        void setValue2  (double value)  {BCValue2 = value;};
        //GETTERS
        int     getType()   {return BCType;};
        double  getValue1() {return BCValue1;};
        double  getValue2() {return BCValue2;};
};


/***************************************************************************************************
 * Input parameters are stored in this class. readSettingsFile() method does the job: It opens the
 * settings.in file, reads the parameters and assigns them to the private variables.
***************************************************************************************************/
class setting
{
    private:
        // VARIABLES
        int     ne;         // number of elements
        int     nn;         // number of nodes
        string  title;      // title of the document
        string  minfFile;   // information file name
        string  mxyzFile;   // node coordinates file name
        string  mienFile;   // connctivity file name
        string  mrngFile;   // boundary info file name
        string  dataFile;   // data file name
        double  initT;      // initial value of the temperature
        double  D;          // Diffusion coefficient
        double  source;     // Heat source term/ force, we call it.
        int     nFG;        // Number of face groups
        int     nIter;      // number of maximum time steps
        double  dt;         // time step size
        bool    restart;    // is this a restart?
        bndc*   BC;         // array for boundary condition types and values

        // METHODS
        void readSettingsFile();    // Reads settings.in file
        void readMinfFile();        // Reads minf file

    protected:

    public:
        // DEFAULT CONSTRUCTOR //
        setting();      //see setting.cpp

        // GETTERS //
        int     getNe()         {return ne;};
        int     getNn()         {return nn;};
        string  getTitle()      {return title;};
        string  getMinfFile()   {return minfFile;};
        string  getMxyzFile()   {return mxyzFile;};
        string  getMienFile()   {return mienFile;};
        string  getMrngFile()   {return mrngFile;};
        string  getDataFile()   {return dataFile;};
        double  getInitT()      {return initT;};
        double  getD()          {return D;};
        double  getSource()     {return source;};
        int     getNFG()        {return nFG;};
        bndc*   getBC(int i)    {return &BC[i];};
        int     getNIter()      {return nIter;};
        double  getDt()         {return dt;};
        bool    getRestart()    {return restart;};

        //INTERFACE METHOD
        void readSettings();
};

#endif /* SETTING_H_ */
