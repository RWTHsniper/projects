/*
 * vtks.cpp
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#include "vtks.h"


vtks::vtks(size_t* argN,double* argdx) {
	// TODO Auto-generated constructor stub
	for(size_t i=0;i<3;i++)	{
		N[i]=argN[i];
		dx[i]=argdx[i];
	}
}

void vtks::init_r_grid_vtkfile(string vtkfilename){
	std::cout << "writing vtk file as 2D rectilinear grid "<< vtkfilename << std::endl;
	std::string	tmpname = vtkfilename + ".vtk";
	outfile.open(tmpname.c_str(), std::fstream::out | std::fstream::trunc);
	if (outfile.is_open()!=1){
		cout << "VTK file could not be created. aborting..."<<endl;
		abort();
	}
	outfile << "# vtk DataFile Version 3.0" <<std::endl<<"vtk output"<<std::endl<<"ASCII"<<std::endl<<"DATASET RECTILINEAR_GRID" << std::endl;
	outfile << "DIMENSIONS " << N[0]<<" "<<N[1]<<" "<<N[2]<<std::endl;
	outfile << "X_COORDINATES "<<N[0]<< " float"<<std::endl;

	cout << "Ns "<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
	cout << "dx,dy,dz "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;

	for (size_t i=0; i< N[0]; i++){outfile << double(i)*dx[0]<<" ";				}
		outfile << std::endl<<"Y_COORDINATES "<<N[1]<< " float"<<std::endl;
				for (size_t i=0; i< N[1]; i++){
						outfile << double(i)*dx[1]<<" ";
					}
			outfile << std::endl<< "Z_COORDINATES "<<N[2]<< " float"<<std::endl;
				for (size_t i=0; i< N[2]; i++) {
					outfile << double(i)*dx[2] <<" ";
			}

	outfile << endl;
}

void vtks::init_r_grid_point_data(){
	outfile << "POINT_DATA "<<N[0]*N[1]*N[2]<<endl;
}


void vtks::close(){
	cout << "closing outfile"<<endl;
    outfile.close();
}

/*
Have a look on header file.
template<typename T>
void vtks::write_point_data<T>(T* point_data, std::string data_name,unsigned int dim_data){
if (dim_data == 1){
// Scalar point data
	outfile << "POINT_DATA "<<N[0]*N[1]*N[2]<<endl;
	outfile << "SCALARS "<< data_name <<" ";
	outfile << "int" << endl;

	if ((typeid(T)==typeid(unsigned int))||(typeid(T)==typeid(int))){
//	Integer data write
		outfile << "int" << endl;
	}
	else{
	//	float data write
		outfile << "float" << endl;
	}

	outfile << "LOOKUP_TABLE default"<<endl;

	helpers::indexing ind(N,3);
	for (size_t i=0; i< N[0]; i++)
		for (size_t j=0; j< N[1]; j++)
			for (size_t k=0; k< N[2]; k++) {
				outfile << point_data[ind.g(i,j,k)]<<std::endl;
			}
}
}
*/

vtks::~vtks() {
	// TODO Auto-generated destructor stub
}

