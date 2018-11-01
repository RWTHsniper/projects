#include <string.h>
#include <stdlib.h>
#include <stdio.h>
# include <fcntl.h>
#include "interfaces.h"

// Use this routine in order to swap bytes
void swapbytes(char* array, int nelem, int elemsize) {
  register int sizet, sizem, i, j;
  char* bytea;
  sizet = elemsize;
  sizem = sizet - 1;
  bytea = malloc(sizet);
  for(i = 0; i < nelem; ++i) {
    memcpy((void*) bytea, (void*) (array + i*sizet), sizet);
    for(j = 0; j < sizet; ++j) {
      array[i*sizet + j] = bytea[sizem - j];
    }
  }
  free(bytea);
}

// Use this routine in order to read a file containing coordinates (mxyz)
void read_nodes(
    double* mxyz,     // Output: An array containing the coordinates
    unsigned nn,      // Input: The number of nodes to read
    unsigned nsd,     // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
    ){
  // YOUR CODE STARTS HERE
int fd;
  // I need to change it
fd = open( file_name , O_RDONLY );
read( fd , mxyz , 8* nsd * nn ); 
swapbytes((char*)mxyz , nsd * nn , 8) ; 
close( fd );
}

// Use this routine in order to read a file containing data (pres)
void read_data(
    double* data,     // Output: An array containing the data
    unsigned nn,      // Input: The number of nodes to read
    unsigned ndf,     // Input: The number of degrees of freedom
    char* file_name   // Input: The path to the data file
    ){
  // YOUR CODE STARTS HERE
int fd;
  // I need to change it
fd = open( file_name , O_RDONLY );
read( fd , data , 8* ndf * nn ); 
swapbytes((char*)data , ndf * nn , 8) ; 
close( fd );
}

// Use this routine in order to read the connectivity file (mien)
void read_connectivity(
    int* mien,        // Output: An array containing the connectivity information
    unsigned ne,      // Input: The number of elements
    unsigned nen,     // Input: The number of nodes per element
    char* file_name   // Input: The path to the connectivity file
    ){
  // YOUR CODE STARTS HERE
int fd;
  // I need to change it
fd = open( file_name , O_RDONLY );
read( fd , mien , 4* ne * nen ); 
swapbytes((char*)mien , ne * nen , 4) ; 
close( fd );
}

// Use this routine in order to write the projected data (pres)
void write_data(
    double* data,     // Input: Data to be written
    unsigned nn,      // Input: The number of nodes to write
    unsigned ndf,     // Input: The number of degrees of freedom
    char* file_name   // Input: The path to the data file
    ){
  // YOUR CODE STARTS HERE

FILE *fd;
fd=fopen(file_name, "wb");
swapbytes((char *)data, nn * ndf, 8);
fwrite(data, sizeof(double), nn * ndf,fd);
close(fd);
}

