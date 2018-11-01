/*
 * fini_diff.h
 *
 *  Created on: Nov 10, 2017
 *      Author: jaeyong
 */

#ifndef FINI_DIFF_H_
#define FINI_DIFF_H_

#include <iostream>     // std::cout
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Sparse>

//#include "sparmat.h"

namespace FD{
size_t	global_indx(size_t i,size_t j,size_t k,size_t* N);
//void	build_CD_matrix(spar::spar_mat& fdmat,size_t* N,double* dx,double alph, unsigned int model);
void	build_CD_matrix(Eigen::SparseMatrix<double,Eigen::RowMajor>& argMat,size_t* N,double* dx,double alph, unsigned int model);
}


#endif /* FINI_DIFF_H_ */
