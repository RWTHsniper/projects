#include "utility.h"

namespace utility
{
    /***********************************************************************************************
     * void swapBytes(char* array, int nelem, int elsize)
     ***********************************************************************************************
     * Byte swapping for endianness.
     **********************************************************************************************/
    void swapBytes (char *array, int nelem, int elsize)
    {
    register int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
        sizem = sizet - 1;
    bytea = new char [sizet];
        byteb = new char [sizet];
    for (i = 0; i < nelem; i++)
        {
        memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++) 
                byteb[j] = bytea[sizem - j];
        memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }
    free(bytea); 
        free(byteb);

    return;
    }

}

