/*
 * grid2d.h
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#ifndef GRID2D_H_
#define GRID2D_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>      // std::ifstream
#include <numeric>      // std::inner_produc
#include <math.h>
#include <algorithm>
#include <set>	//	used for sorting..

//	My hearder
#include "indexing.h"

using namespace std;

class grid_2d {
private:
//	unsigned int N[3];
	double L[3];
//	vector<double> xyz;
	vector<unsigned int> mat_indx;

public:
	grid_2d(size_t Nx, size_t Ny, size_t Nz,double Lx,double Ly,double Lz);
	virtual ~grid_2d();

//	public variables
	unsigned int num_materials;
	size_t N[3];
	double	dx[3];
//	std::vector<double> dx(3);

//	Public methods
	void make_homo_grid(unsigned int mat_indx);
	void read_grid(std::string gridfilename, char type);
//	void summary();
	unsigned int* get_mat_indx(){return &mat_indx[0];};
	void summary();
//	double* get_xyz(){return &xyz[0];};

};

#endif /* OTHER_SOURCE_CODES_GRID2D_H_ */
