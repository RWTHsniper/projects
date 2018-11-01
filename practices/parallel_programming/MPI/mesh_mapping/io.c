#include <string.h> // memcpy
#include <stdlib.h> // malloc
#include <stdio.h>  // printf

#include <unistd.h> // lseek,read,write
#include <fcntl.h> // open etc.

#include <err.h>

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

// Use this routine in order to read a file containing double data
void read_data(
    double* data,     // Output: Array containing the data
    size_t nn,        // Input: Number of nodes/elements to read
    off_t offset,     // Input: Number of nodes/elements to offset the reading
    unsigned entries, // Input: Number of data entries (DOFs or coordinates) per node/element
    char* file_name   // Input: Path to the data file
    ){
  int fd;
  fd = open(file_name, O_RDONLY);
  if (fd == -1) {
    err(1,"Opening %s for reading", file_name);
  }

// code here
  lseek(fd, offset*entries*sizeof(double), SEEK_SET);

  size_t bytestoread = nn*entries*sizeof(double);// CODE HERE: Number of bytes to read
  size_t bytesread = read(fd, data, bytestoread);
  if (bytesread == -1) {
    err(1, "%s: Reading data failed", file_name);
  } else if (bytesread != bytestoread) {
    printf("%s: Only %zd bytes of %zu bytes read!\n", file_name, bytesread, bytestoread);
    exit(1);
  }

  close(fd);

  swapbytes((char*)data, nn*entries, sizeof(double));

  // DO NOT CHANGE THE FOLLOWING LINE
  check_input(data,nn,offset,entries,file_name);
}


// Use this routine in order to write the projected data
void write_data(
    double* data,     // Input: Data to be written
    size_t nn,        // Input: Number of nodes to write
    off_t offset,     // Input: Number of elements to offset
    unsigned entries, // Input: Number of degrees of freedom
    char* file_name   // Input: Path to the data file
    ){
  // DO NOT CHANGE THE FOLLOWING LINE
  check_output(data,nn,offset,entries,file_name);

  swapbytes((char*)data, nn*entries, sizeof(double));
  int fd;
  fd = open(file_name, O_RDWR | O_CREAT, 0666);
  if (fd == -1) {
    err(1,"Opening %s for writing", file_name);
  }

  lseek(fd, offset*entries*sizeof(double), SEEK_SET);

  size_t bytestowrite = nn*entries*sizeof(double);// CODE HERE: Number of bytes to read
  size_t byteswritten = write(fd, data, bytestowrite);
  if (byteswritten == -1) {
    err(1, "%s: Writing data failed", file_name);
  } else if (byteswritten != bytestowrite) {
    printf("%s: Only %zd bytes of %zu bytes written!\n", file_name, byteswritten, bytestowrite);
    exit(1);
  }

  close(fd);
}


