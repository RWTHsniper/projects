/*
 * grid2d.cpp
 *
 *  Created on: Nov 1, 2017
 *      Author: jaeyong
 */

#include "grid2d.h"
#define PI 3.14159265
using namespace std;

//	It is like consideering dimension as 3D.

grid_2d::grid_2d(size_t Nx, size_t Ny, size_t Nz,double Lx,double Ly,double Lz) {
	// TODO Auto-generated constructor stub
N[0]=Nx;N[1]=Ny;N[2]=Nz;L[0]=Lx;L[1]=Ly;L[2]=Lz;
mat_indx.resize(Nx*Ny*Nz,0);
num_materials=1;	// By default it is 1. It can be changed when composite RVE is read.
//xyz.resize(Nx*Ny*3,0);	// xyz[Nx][Ny][2]
dx[0] = Lx/double(Nx);
dx[1] = Ly/double(Ny);
dx[2] = Lz/double(Nz);
/*
for (size_t i=0;i<Nx;i++)
	for(size_t j=0;j<Ny;j++)
		for (size_t k=0;k<Nz;k++){
		xyz[i*Ny*3+j*3+0] = dx[0]*double(i);
		xyz[i*Ny*3+j*3+1] = dx[1]*double(j);
		xyz[i*Ny*3+j*3+2] = dx[2]*double(k);
	}
	*/
}

void grid_2d::summary(){
for (size_t i=0;i<40;i++) cout << "*";
cout << endl<< "Summary of grid information"<<endl;
cout<<"Size of grid space: ["<<N[0]<<" "<<N[1]<<" "<<N[2]<<"]"<<endl;
cout<<"Length of grid space: ["<<L[0]<<" "<<L[1]<<" "<<L[2]<<"]"<<endl;
cout<<"dx of grid space: ["<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<"]"<<endl;
cout<<"Number of material types: ["<<num_materials<<"]"<<endl;
if (num_materials==1)	cout<<"Homogeneous RVE"<<endl;
else	cout<<"Composite RVE"<<endl;
for (size_t i=0;i<40;i++) cout << "*";
cout << endl;

}

void grid_2d::make_homo_grid(unsigned int homo_indx){
	helpers::indexing ind(N,3);
	for (size_t i=0;i<N[0];i++)
		for (size_t j=0;j<N[1];j++)
			for (size_t k=0;k<N[2];k++){
				mat_indx[ind.g(i,j,k)] = homo_indx;
		}

}


void grid_2d::read_grid(std::string gridfilename, char type){

	//	Prestep:	Start counting number of lines in a file.
	std::ifstream dummyfile;
	dummyfile.open(gridfilename.c_str());
	size_t num_lines=0;
	string dummyline;
    while (std::getline(dummyfile, dummyline)) ++num_lines;
	dummyfile.close();

	if	(N[0]*N[1]*N[2] != num_lines){
		cout << "Size of grid and file mismatch!!" << endl;
		cout << "number of lines in grid file "<<num_lines<<endl;
		cout<<"Size of grid space: ["<<N[0]<<" "<<N[1]<<" "<<N[2]<<"] = "<< (N[0]*N[1]*N[2]) <<endl;
		exit(-1);
	}

	std::ifstream gridfile;
	gridfile.open(gridfilename.c_str());
	std::cout<<"start reading grid file "<<gridfilename<<std::endl;


	helpers::indexing ind(N,3);
	if (gridfile.is_open()){
		if (type == 'C'){
			for (size_t i=0;i<N[0];i++)
				for (size_t j=0;j<N[1];j++)
					for (size_t k=0;k<N[2];k++){
						gridfile >> mat_indx[ind.g(i,j,k)];
	//					printf("mat_indx, %d \n",mat_indx[i][j][k]);
				}
		}
		else if (type == 'F'){
			for (size_t k=0;k<N[2];k++)
				for (size_t j=0;j<N[1];j++)
					for (size_t i=0;i<N[0];i++){
						gridfile >> mat_indx[ind.g(i,j,k)];
	//					printf("mat_indx, %d \n",mat_indx[i][j][k]);
				}
		}

// using default comparison:
/*
		cout << "print mat_indx array"<<endl;
		for (size_t i=0;i<(N[0]*N[1]*N[2]);i++) std::cout << mat_indx[i];
*/
		/*
		set<unsigned int> sa(&mat_indx[0],&mat_indx[N[0]*N[1]*N[2]]);
//		  for (size_t i=0;i<sa.size();i++) std::cout << sa[i];
		  cout<<endl;
		  std::cout << "size of unique material properties " << sa.size() << std::endl;
		  num_materials=sa.size();
		  sa.clear();
		  */

	  std::vector<unsigned int> myvect(&mat_indx[0],&mat_indx[N[0]*N[1]*N[2]]);
	  sort( myvect.begin(), myvect.end() );
	  myvect.erase( unique( myvect.begin(), myvect.end() ), myvect.end() );
	  std::cout << "size of unique material properties " << myvect.size() << std::endl;
	  std::cout << "Contents of material indexes " << std::endl;
	  for (size_t i=0;i<myvect.size();i++) std::cout << myvect[i]<< " ";
	  cout<<endl;
	  num_materials=myvect.size();
	  myvect.clear();
	}
	else{
		std::cout<<"gridfile is not opened"<<std::endl; abort();
	}
	gridfile.close();
}


grid_2d::~grid_2d() {
	// TODO Auto-generated destructor stub
}

