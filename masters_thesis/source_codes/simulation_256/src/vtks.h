/*
 * vtks.h
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#ifndef OTHER_SOURCE_CODES_VTKS_H_
#define OTHER_SOURCE_CODES_VTKS_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <typeinfo>
#include "indexing.h"


using namespace std;


class vtks {
private:
	size_t N[3];
	double	dx[3];
	std::ofstream outfile;
public:
	vtks(size_t* argN,double* argdx);
	virtual ~vtks();
//	Public methods
	void init_r_grid_vtkfile(std::string vtkfilename);
	void close();

//	members using template should be defined in hearder file.
	/*	write_point_data
	 * If dim_data == 1 -> scalar
	 * if dim_data > 1	-> vector
	 */
/*
 * init_r_grid_point_data:
 * Initialized rectilinear grid vtk file.
 * This function sets up Dimension and point coordinates of vtk file.
 */
	void init_r_grid_point_data();
/*
 * This member functions to write point data on initialized vtk file.
 * init_r_grid_vtkfile should be called first.
 */
	template<typename T>
	void write_scalar_point_data(T* point_data, std::string data_name,unsigned int dim_data){
	// Scalar point data
		outfile << "SCALARS "<< data_name <<" ";
		if ((typeid(T)==typeid(size_t))||(typeid(T)==typeid(unsigned int))||(typeid(T)==typeid(int))){
	//	Integer data write
			outfile << "int "<<dim_data << endl;
		}
		else{
		//	float data write
			outfile << "float "<<dim_data << endl;
		}
		outfile << "LOOKUP_TABLE default"<<endl;
		size_t Nind[] = {N[0],N[1],N[2],dim_data};
		helpers::indexing ind(Nind,4);
/*
		for (size_t i=0; i< N[0]; i++)
			for (size_t j=0; j< N[1]; j++)
				for (size_t k=0; k< N[2]; k++){
					for (size_t l=0;l<dim_data;l++){
					outfile << point_data[ind.g(i,j,k,l)]<<" ";
					}
					outfile << endl;
				}
*/
		for (size_t k=0; k< N[2]; k++)
			for (size_t j=0; j< N[1]; j++)
				for (size_t i=0; i< N[0]; i++){
					for (size_t l=0;l<dim_data;l++){
					outfile << point_data[ind.g(i,j,k,l)]<<" ";
					}
					outfile << endl;
				}

	}

	template<typename T>
	void write_vector_point_data(T* point_data, std::string data_name){
	// Vector point data
		outfile << "VECTORS "<< data_name <<" ";
		if ((typeid(T)==typeid(unsigned int))||(typeid(T)==typeid(int))){
	//	Integer data write
			outfile << "int "<< endl;
		}
		else{
		//	float data write
			outfile << "float "<< endl;
		}
		unsigned int dim_data = 3;
		size_t Nind[] = {N[0],N[1],N[2],dim_data};
		helpers::indexing ind(Nind,4);
	/*
		for (size_t i=0; i< N[0]; i++)
			for (size_t j=0; j< N[1]; j++)
				for (size_t k=0; k< N[2]; k++){
					for (size_t l=0;l<dim_data;l++){
					outfile << point_data[ind.g(i,j,k,l)]<<" ";
					}
					outfile << endl;
				}
				*/
		for (size_t k=0; k< N[2]; k++)
			for (size_t j=0; j< N[1]; j++)
				for (size_t i=0; i< N[0]; i++){
					for (size_t l=0;l<dim_data;l++){
					outfile << point_data[ind.g(i,j,k,l)]<<" ";
					}
					outfile << endl;
				}
	}


};

#endif /* OTHER_SOURCE_CODES_VTKS_H_ */
