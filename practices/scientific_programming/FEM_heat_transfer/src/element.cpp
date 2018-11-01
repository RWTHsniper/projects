#include "element.h"
#include "utility.h"
/***************************************************************************************************
 * void element::createElements()
 * This routine calls private methods to read mien and mrng file.
 **************************************************************************************************/
void element::createElements(setting* settings, node* nodes)
{
    // Initialize all private variables
    for (int i = 0; i < ne*nen; i++)
    {
	conn[i] = 0;
	M[i]    = 0.0;
	F[i]    = 0.0;
    }
    for (int i = 0; i < ne*nef;     i++)    FG[i]  = 0;
    for (int i = 0; i < ne*nen*nen; i++)    K[i]   = 0.0;
    // Start element operations
    readMien(settings);
    readMrng(settings);
    calculateJacobian(nodes);
    calculateElementMatrices(settings, nodes);
    applyBoundaryConditions(settings, nodes);
    return;
}
/***************************************************************************************************
 * void element::readMien()
 * This routine reads mien file.
 **************************************************************************************************/
void element::readMien(setting* settings)
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files
    dummy = settings->getMienFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
	cout << "Unable to open file : " << dummy << endl;
	exit(0);
    }
    readStream = new char [nen*sizeof(int)];
    file.seekg (0, ios::beg);
    for(size_t i=0; i<ne; i++)
    {
	file.read (readStream, nen*sizeof(int));
	utility::swapBytes(readStream, nen, sizeof(int));
	for(int j=0; j<nen; j++)
	    conn[i*nen+j] = *((int*)readStream+j)-1;
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();
    return;
}
/***************************************************************************************************
 * void element::readMrng()
 * This routine reads mrng file.
 **************************************************************************************************/
void element::readMrng(setting* settings)
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files
    dummy = settings->getMrngFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
	cout << "Unable to open file : " << dummy << endl;
	exit(0);
    }
    readStream = new char [nef*sizeof(int)];
    file.seekg (0, ios::beg);
    for(size_t i=0; i<ne; i++)
    {
	file.read (readStream, nef*sizeof(int));
	utility::swapBytes(readStream, nef, sizeof(int));
	for(ptrdiff_t j=0; j<nef; j++)
	{
	    FG[i*nef+j] = *((int*)readStream+j); 
	    //cout << FG[i*nef+j] << endl;
	}
    }
    cout << "> File read is completed: " << dummy << endl;
    file.close();
    return;
}
/***************************************************************************************************
 * void element::calculateJacobian()
 * This routine calculates the Jacobian determinant and derivatives for each GQP.	det(J) / derivatives for each GQP
 **************************************************************************************************/
void element::calculateJacobian(node* nodes)
{
    // look for (2) in Ex2
//    double J[nsd*nsd];             // J[2*2]
    double J[nsd][nsd];             // J[2*2]	
///    double A[nsd*nen], B[nen*nsd];  // dSdxi & dS/deta, node array  
    double A[nsd][nen], B[nen][nsd];  // dSdxi & dS/deta, node array  
//    double inv_J[nsd*nsd];
    double inv_J[nsd][nsd];	
    for (size_t ie=0; ie<ne;ie++)  // index of element in global
	{
	    for (size_t ien=0; ien<nen;ien++) // index of element node
	    {   // B matrix is element's node array
		B[ien][0]   = nodes->getX(conn[nen*ie+ien]);
		B[ien][1] = nodes->getY(conn[nen*ie+ien]); 
		}		
	for (size_t iGQP=0; iGQP<nGQP;iGQP++) // index og GQP
	{
	    for (size_t ien=0; ien<nen;ien++) // index of element node
	    {   // Assemble A in column, B in row
		A[0][ien]       = GaussQ::getDSdKsi(iGQP,ien); 
		A[1][ien]  = GaussQ::getDSdEta(iGQP,ien);
	    }
	    // Assemble Jacobian matrix
	    for (size_t k=0; k<nsd;k++)
	    {
		for (size_t l=0;l<nsd;l++)
		{   // J[k][l] = A[k-th row]*B[l-th column]
			J[k][l]=A[k][0]*B[0][l]+A[k][1]*B[1][l]+A[k][2]*B[2][l];
		    //cout << "Jacobian" << (nsd*k+l) << J[nsd*k+l] << endl;
		}
	    }
	    double detJ=J[0][0]*J[1][1]-J[0][1]*J[1][0];
	    GaussQ::setDetJ(ie,iGQP,detJ);		
	    inv_J[0][0]=J[1][1]/detJ;
	    inv_J[0][1]=-J[0][1]/detJ;
	    inv_J[1][0]=-J[1][0]/detJ;
	    inv_J[1][1]=J[0][0]/detJ;
	    for(size_t ien=0; ien<nen; ien++)
	    {
		GaussQ::setDSdx(ie,iGQP,ien,inv_J[0][0]*GaussQ::getDSdKsi(iGQP,ien)+inv_J[0][1]*GaussQ::getDSdEta(iGQP,ien));
		GaussQ::setDSdy(ie,iGQP,ien,inv_J[1][0]*GaussQ::getDSdKsi(iGQP,ien)+inv_J[1][1]*GaussQ::getDSdEta(iGQP,ien));
	    }
	}
	}
    return;
}
/***************************************************************************************************
 * void element::calculateElementMatrices()
 * This routine calculates element mass, stiffness and forcing matrices
 * as well as nodal mass.
 * Global matrix size
 conn[ne*nen]
 M[ne*nen]
 F[ne*nen]
 FG[ne*nef]
 K[ne*nen*nen]
 **************************************************************************************************/
void element::calculateElementMatrices(setting* settings, node* nodes)
{    
    for (size_t ie=0; ie<ne;ie++)
    {// Setting element mass terms
	double M_elm[nen][nen];
	for (size_t ien=0; ien<nen;ien++)
	{
	    for (size_t jen=0; jen<nen;jen++)
	    {
		M_elm[ien][jen]=0;
	    }
	}
	double M_e=0;
	double M_tr=0;
	for (size_t ien=0; ien<nen; ien++)
	{   // calculate element mass matrix
	    for (size_t jen=0;jen<nen;jen++)
	    {
		for (size_t iGQP=0;iGQP<nGQP;iGQP++)
		{
		    double Si=GaussQ::getS(iGQP,ien);
		    double Sj=GaussQ::getS(iGQP,jen);
		    double weight = GaussQ::getWeight(iGQP);
		    double detJ = GaussQ::getDetJ(ie,iGQP);
			M_elm[ien][jen] += Si*Sj*weight*detJ;
		}
	    }
	}
	for (size_t ien=0;ien<nen;ien++)
	{// calculate M_e and M_tr to get Lumped mass system
	    M_tr += M_elm[ien][ien];
	    for (size_t jen=0;jen<nen;jen++) M_e += M_elm[ien][jen];
	}
	for (size_t ien=0;ien<nen;ien++)
	{// calculate M_elm_lump
	    double M_lump=M_elm[ien][ien]*M_e/M_tr;
	    setM(ie,ien,M_lump);
	}
    //}
    // calculate Nodal mass as well
    // nodes->addMass;
    //for (size_t ie=0; ie<ne;ie++)
	//{	
		for (size_t ien=0; ien<nen; ien++)
	{
	    size_t i=conn[ie*nen+ien];
	    double mass = M[ie*nen+ien];
	    nodes->addMass(i, mass);
	}
	//}

    // Set up element K_ij
    double k=settings->getD();      
    //for (size_t ie=0; ie<ne;ie++)
    //{// Loop over each elem
	for (size_t ien=0;ien<nen;ien++)
	{
	    for (size_t jen=0;jen<nen;jen++)
	    {
		double K_temp=0;
		for (size_t iGQP=0;iGQP<nGQP;iGQP++)
		{// Integrate using GQP
		    double dSi_x=GaussQ::getDSdx(ie,iGQP,ien);
		    double dSi_y=GaussQ::getDSdy(ie,iGQP,ien);   
		    double dSj_x=GaussQ::getDSdx(ie,iGQP,jen);
		    double dSj_y=GaussQ::getDSdy(ie,iGQP,jen);    
		    double detJ=GaussQ::getDetJ(ie,iGQP); 
		    double weight=GaussQ::getWeight(iGQP);
		    K_temp += detJ*((dSi_x*dSj_x)+(dSi_y*dSj_y))*weight;
		}	// we don't multiply k
		setK(ie,ien,jen, K_temp);   // Set K matrix on element level
	    }
	}
    }
    //Forcing matrices
    // Setting element force terms
    for (size_t ie = 0; ie < ne;ie++){
		// Loop over each elem
	for (size_t ien = 0; ien < nen; ien++){   
	    double temp = 0;
	    for (size_t iGQP = 0; iGQP < nGQP; iGQP++){
		// calculate element mass matrix
		double weight=GaussQ::getWeight(iGQP);
		double Si = GaussQ::getS(iGQP, ien);
		double Source = settings->getSource();
		double detJ = GaussQ::getDetJ(ie,iGQP);
		temp += weight*Si*Source*detJ; 
	    }
	    setF(ie, ien, temp);
	}
    }
    return;
}
/***************************************************************************************************
 * void element::applyBoundaryConditions()
 * Applies Drichlet and Neumann boundary conditions.
 **************************************************************************************************/
void element::applyBoundaryConditions(setting* settings, node* nodes)
{
    // Read FG and check BC type
    // Specify BC type on two nodes
    // Calculate Neumann or Dirichlet BC
    for (size_t ie=0; ie<ne;ie++)
    {
	for (size_t iefTri=0;iefTri< nefTri; iefTri++)
	{
	    int f_Dc[2] = {0,0};  // Flags which already has Dirichlet BC
	    int face_code=FG[ie*nef+iefTri];
	    if (face_code>0)
	    {// We found BC. iefTri: face_index
		int node[2]= {edgeNodesTri[iefTri][0], edgeNodesTri[iefTri][1]};
		int g_node[2]= {conn[ie*nen+node[0]], conn[ie*nen+node[1]]};
		// set BCtype
		bndc* temp_bndc=settings->getBC(face_code);	
		int typ = temp_bndc->getType();	// 1: Dirichlet, 2: Neumann
		if ((nodes->getBCtype(g_node[0]) > 0) || (nodes->getBCtype(g_node[1]) > 0))   
		{
//		    cout << "Multiple boundary condition detected" << endl;
		    if (typ == 2)  
		    {   // Neumann BC will be applied only on a node without Dirichlet BC.
			if (nodes->getBCtype(g_node[0]) == 1)
			{
			    f_Dc[0] = 1;
			}
			if (nodes->getBCtype(g_node[1]) == 1)
			{
			    f_Dc[1] = 1;
			}
		    }
		}
		if (typ == 1)
		{
		    // Apply dirichlet BC
		    double bc_v = temp_bndc->getValue1();
		    nodes->setBCtype (g_node[0], typ); nodes->setBCtype (g_node[1], typ);
		    nodes->setT(g_node[0], bc_v); nodes->setT(g_node[1], bc_v);
		}
		else if (typ == 2)   // Neumann BC
		{
		    if ((f_Dc[0] == 1) && (f_Dc[1] == 1))
		    { 
			continue;   // If both nodes have Dirichlet BC already, skip it.
		    }
		    else
		    { 
			double bc_v = temp_bndc->getValue1();
			double x[2] ={nodes->getX(g_node[0]), nodes->getX(g_node[1])};
			double y[2] ={nodes->getY(g_node[0]), nodes->getY(g_node[1])};
			double dx = x[1]-x[0]; double dy = y[1]-y[0]; 
			double L= sqrt(dx*dx+dy*dy);    // length of an face
			double q_in = bc_v;				
			double B_value = q_in*L/2;
			// Apply Neumann BC only if a node does not have Dirichlet BC
			if (f_Dc[0] != 1) 
			{
			    nodes->setBCtype (g_node[0], typ);
			    addF (ie, node[0], B_value); 
			}				
			if (f_Dc[1] != 1)
			{
			    nodes->setBCtype (g_node[1], typ);
			    addF (ie, node[1], B_value);
			}
		    }
		}
	    }
	}
    }
    return;
}