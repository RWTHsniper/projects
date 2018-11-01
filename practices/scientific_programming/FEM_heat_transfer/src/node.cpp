#include "node.h"
#include "utility.h"

/***************************************************************************************************
 * void node::createNodes()
 * This routine calls private methods to read mxyz and data file.
 **************************************************************************************************/
void node::createNodes(setting* settings)
{
    // Initialize the private variables with zero.
    for (int in = 0; in < nn; in++)
    {
        BCtype[in] = 0;
        x[in]      = 0.0;
        y[in]      = 0.0;
        T[in]      = 0.0;
        mass[in]   = 0.0;
    }

    readMxyz(settings);
    readData(settings);
    
    return;
}


/***************************************************************************************************
 * void node::readMxyz()
 ***************************************************************************************************
 * READ THE MXYZ FILE
 * mxyz file contains the node coordinates
 **************************************************************************************************/
void node::readMxyz(setting* settings)
{
    int         nn;             // nuber of nodes
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files

    nn    = settings->getNn();
    dummy = settings->getMxyzFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    readStream = new char [nsd*sizeof(double)];
    file.seekg (0, ios::beg);
    for(int i=0; i<nn; i++)
    {
        file.read(readStream, nsd*sizeof(double));
        utility::swapBytes(readStream, nsd, sizeof(double));
        x[i] = *((double*)readStream);
        y[i] = *((double*)readStream+1);
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();

    return;
}


/***************************************************************************************************
 * void node::readData()
 ***************************************************************************************************
 * READ THE DATA FILE
 * data file contains the initial value of the field.
 **************************************************************************************************/
void node::readData(setting* settings)
{
    int      nn;             // number of nodes to be determined from minf file
    bool     restart;        // Is this a restart?
    ifstream file;           // file name obj
    string   dummy;          // dummy string to hold names
    char*    readStream;     // temperory var used for strings read from files
    double   dummyDouble;    // temperory var used for double values read from files

    nn      = settings->getNn();
    restart = settings->getRestart();

    if (restart)
    {
        dummy = settings->getDataFile();
        file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
        if (file.is_open()==false)
        {
            cout << "Unable to open file : " << dummy << endl;
            exit(0);
        }
        readStream = new char [sizeof(double)];
        file.seekg (0, ios::beg);
        for(int i=0; i<nn; i++)
        {
            file.read (readStream, sizeof(double));
            utility::swapBytes(readStream, 1, sizeof(double));
            T[i] = *((double*)readStream);
        }
        cout << "> File read is completed: " << dummy << endl;
        file.close();
    }
    else
    {
        dummyDouble = settings->getInitT();
        for(int i=0; i<nn; i++)
            T[i] = dummyDouble;
        cout << "> Data initialization is completed." << endl;
    }

    return;
}
