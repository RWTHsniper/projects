#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 09:48:57 2018

@author: jungjaeyong

Momentum balance equation problem



"""

"""
Neo_Hook

variables
    F: deformation gradient
    C: C=F'*F right cauchy-green tensor
    E: E=0.5(C-I) green-lagrangian tensor

"""

"""
    .bc = {'d':,'n':}
    'd': dirichlet bc. (displacement) made of [nodeindex, dof, value]
    'n': neumann bc. (traction)
"""

"""
MBE
    outp_lst: list of output field quantities
    activ_dof: map tot_dof to activ_dof
    inv_glob_dof: map acti_dof to tot_dof

"""

import math
import numpy as np
Dtyp = np.float64
Ityp = np.int64
import GQ
import MBE_mats
import helpers

#import mesh

    
        
class BC:
    def __init__(self):
        self.bcs = {'d':[],'p':{},'n':[]}

#   add dirichlet_bc
#   index: index of node. (it can be a list of indexes)
#   direction: o(x), 1(y), 2(z)
    def add_dirichlet_bc(self,indx,direction,val):
        if (type(direction)==str):
            if (direction=="x"):
                direction = 0
            elif (direction=="y"):
                direction = 1
            elif (direction=="z"):
                direction = 2
        if (type(indx)==list or type(indx)==np.array or type(indx)==np.ndarray):
            for i in range(len(indx)):
                self.bcs['d'].append([indx[i],(direction),Dtyp(val)])
        else:
            self.bcs['d'].append([indx,(direction),Dtyp(val)])

    def replace_dirichlet_bc(self,indx,direction,val):
        if (type(direction)==str):
            if (direction=="x"):
                direction = 0
            elif (direction=="y"):
                direction = 1
            elif (direction=="z"):
                direction = 2
        if (type(indx)==list or type(indx)==np.array or type(indx)==np.ndarray):
            for i in range(len(indx)):
                flag = 0
                for b_indx,d_lst in enumerate(self.bcs['d']):
                    if (indx[i]==d_lst[0] and direction==d_lst[1]):
                        d_lst[2] = Dtyp(val)
                        flag = 1
                        break
                if (flag==0): self.bcs['d'].append([indx[i],(direction),Dtyp(val)])
        else:
            for b_indx,d_lst in enumerate(self.bcs['d']):
                flag = 0
                if (indx==d_lst[0] and direction==d_lst[1]):
                    d_lst[2] = Dtyp(val)
                    flag = 1
                    break
            if (flag==0): self.bcs['d'].append([indx,(direction),Dtyp(val)])
        

#   add periodic boundary condition
#   node_id, node_id_periodic, displacement
#   e.g.) u_10 = u_1 + [10.0,0.0,0.0]. This is represented and saved like..
#   p.update(10,1,np.array([10.0,0.0,0.0],dtype=Dtyp))                
#   There is a rule. Make sure that ind1 > ind2
    def add_perodic_bc(self,ind1,ind2,disp):
        if (ind1 > ind2):
            self.bcs['p'].update({ind1:[ind2,np.array(disp,dtype=Dtyp)]})
        else:
            print ("make sure ind2 smaller than ind1 to assign periodic bc.")
            self.bcs['p'].update({ind2:[ind1,np.array(-disp,dtype=Dtyp)]})

    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")        

"""
MBE: momentum balance equation problem
members of MBE
    mesh: mesh of the MBE problem
    nn: total number of nodes
    bc: boundary condition 
    outp_lst: lists of outputs
    sdim: spatial dimension of the problem. baisically, it is 3-dim
    ndof: total number of dof for a MBE problem
    activ_dof: list of active_dofs
        If a dof is inactive, it has nn
        If a dof is periodic, the corresponding dof is set
    inv_glob_dof: inverse mapping from active value to inactive values
    constitutive_laws
        0: Neo-Hook(purely elastic) = [lambda,mu]
        1: Neo-Hook with isotropic plasticity
        2: Neo-Hook with phenomenological
"""

class MBE:
    def __init__(self,gmesh,**kwargs):
        self.mesh = gmesh
        self.nn = gmesh.nn
        self.bc = BC()
        self.outp_lst=[]
        self.sdim = 3
        self.ndof = self.sdim * self.nn
        self.indx_inactiv = self.ndof
        self.activ_dof = np.zeros(self.ndof,dtype=Ityp)
        self.n_activ_dof = 0
        self.glob_dof = np.zeros(0,dtype=Ityp) # D.O.F. dictionaries of unconstrained D.O.F.s. Mapping glob_dof to solution vector D.O.F.
        self.inv_glob_dof = np.zeros(0,dtype=Ityp) # Mapping D.O.F. of solution vector to the global D.O.F.
        self.solution = np.zeros(0,dtype=Dtyp) # solution vector
#       gq_tet4 = GQ.GQ_tet4() will be implemented in the future        
        self.GQs = {'hex8':GQ.GQ_hex8()}
        self.nel_params = 20 # maximum number of parameters for each contitutive law
        self.nhex8 = self.mesh.elements['hex8'][1]
        self.ntet4 = self.mesh.elements['tet4'][1]
        self.el_params = {'hex8':np.zeros((self.nhex8,self.GQs['hex8'].ngq,self.nel_params),dtype=Dtyp)}
#   field quantities
        self.disp = np.zeros((self.nn,self.sdim),dtype=Dtyp)
#        self.disp = np.ones((self.nn,self.sdim),dtype=Dtyp)
        self.init_field_vari()
        self.mat_type = "dense" #default setting for mat_type is "direct". It can be "sparse"
        self.nlgeom = True
        if (kwargs is not None):
            for key, value in kwargs.items():
                if (key == "mat_type"):
                    self.mat_type = value
                if (key == "nlgeom"):
                    if ((value == True)|(value == False)):
                        self.nlgeom = value
                    else:
                        print("nlgeom should be True or False. By default, nlgeom = True is set.")
        self.init_rKf()
                
#    P: 1st Piola-Kirchhoff stress
#    S: 2nd Piola-Kirchhoff stress
#    F: deformation gradient        
    def init_field_vari(self,*args):
        if (self.mesh.elements['hex8'][1] > 0):
            self.P = {'hex8':np.zeros((self.nhex8,self.GQs['hex8'].ngq,self.sdim,self.sdim),dtype=Dtyp)}
            self.S = {'hex8':np.zeros((self.nhex8,self.GQs['hex8'].ngq,self.sdim,self.sdim),dtype=Dtyp)}
            self.F = {'hex8':np.zeros((self.nhex8,self.GQs['hex8'].ngq,self.sdim,self.sdim),dtype=Dtyp)}
            for elem in args:
                if (elem == "Cauchy"):
                    self.Cauchy = {'hex8':np.zeros((self.nhex8,self.GQs['hex8'].ngq,self.sdim,self.sdim),dtype=Dtyp)}

#   Sort ans check dirichlet B.C.s.
    def sort_check_dirichlet(self):
        self.bc.bcs['d'] = sorted(self.bc.bcs['d'])     
        self.bc.bcs['d'] = helpers.unique(self.bc.bcs['d'])
        try:            
            for indx,elem in enumerate(self.bc.bcs['d']):
                if ((elem[0] == self.bc.bcs['d'][indx+1][0]) and (elem[1] == self.bc.bcs['d'][indx+1][1]) and (elem[2] != self.bc.bcs['d'][indx+1][2])):
                    print("Overlapping DIrichlet B.C.s. with different value!")
                    print(elem)
                    print(self.bc.bcs['d'][indx+1]) 
        except:
            return
                                                     
#   This function should run after BC is defined
    def set_activ_dof(self):
#   Before trying this, sort Dirichlet B.C.s
        self.sort_check_dirichlet()
        for indx in range(self.nn):
            for i in range(self.sdim):
                self.activ_dof[indx*self.sdim+i] = indx*self.sdim + i
#   Dirichlet BC is taken into account        
        for indx,d_bc in enumerate(self.bc.bcs['d']):
            self.activ_dof[d_bc[0]*self.sdim+d_bc[1]] = self.indx_inactiv # inactive!

        for key,value in self.bc.bcs['p'].items():
            self.activ_dof[key*self.sdim:(key+1)*self.sdim] = [value[0]*self.sdim+i for i in range(self.sdim)]
                                    
            
#   Neumann BC will be taken into account someday
            
#   Then, it is time to make inv_glob_dof mapper  
        lst = []    
        for indx,elem in enumerate(self.activ_dof):           
            if (elem == indx): # a dof free of Dirichlet and peroidic B.C.s           
                lst.append(elem)
        self.inv_glob_dof = np.array(lst,dtype=Ityp)
        self.n_activ_dof = self.inv_glob_dof.shape[0]
        self.solution = np.resize(self.solution,self.n_activ_dof)
        self.glob_dof = np.ones(self.ndof,dtype=Ityp) * self.indx_inactiv
# assign mapping of activ_dof to solution vector
        for indx,elem in enumerate(self.inv_glob_dof):
            self.glob_dof[elem] = indx
# mapping of dofs under periodic B.C.
        for key,value in self.bc.bcs['p'].items():
            node = Ityp(key)
            # assign constrained periodic node's dof.
            for i in range(self.sdim):
                self.glob_dof[node*self.sdim+i] = self.glob_dof[value[0]*self.sdim+i]
        self.assign_d_bc()
        self.init_rKf()
        
#   count total number of active dof

    def init_rKf(self):
        if ((self.mat_type == "dense") or (self.mat_type == "d")):
            self.Kg = np.zeros((self.n_activ_dof,self.n_activ_dof),dtype=Dtyp)        
        self.Resi = np.zeros((self.n_activ_dof),dtype=Dtyp)
        self.Fg = np.zeros((self.n_activ_dof),dtype=Dtyp)

#   assign dirichlet bc to displacement array
    def assign_d_bc(self):
        for indx, elem in enumerate(self.bc.bcs['d']):
            self.disp[elem[0]][elem[1]] = elem[2]

    def loop_elem(self):
#   Initialize global stiffness, residual, force tensors        
        self.Kg.fill(0.0)
        self.Fg.fill(0.0)
        self.Resi.fill(0.0)
#   loop over for hex8        
        if (self.nhex8>0):
            self.loop_hex8()
        
            
    def calc_k_elm(self,B_L,dSdE,w,detJ,k_elm):
        k_elm[:,:] += np.matmul(B_L.transpose(),np.matmul(dSdE,B_L)) * w  * detJ

#   Express stress sig using Manel's notation
#   sig = [sig00,sig11,sig22,sqrt(2)*sig01,sqrt(2)*sig12,sqrt(2)*sig20]        
#   fint_i = B^T(i,j)*sig(j)*w
#   fint_i = B(j,i)*sig(j)*w        
    def calc_fint(self,B_L,sig,w,detJ,fint):
        sig2 = MBE_mats.T2_Mandel(sig)
        fint[:] += np.matmul(B_L.transpose(),sig2) * w * detJ

#   assemble global stiffness matrix
#   I think I need time to think about how to apply P.B.C. and calculate global Kg and Fg
    def assemble_Kg_Fg(self,conn,k_el,f_el):
        # dofs in an element
        el_dof = np.zeros((conn.shape[0]*self.sdim),dtype=Ityp)
        for i in range(self.sdim):
            el_dof[i::self.sdim] = (conn*self.sdim + i)
        # loop of an element's indexes
        for indx,idof in enumerate(el_dof):
            mapped_dof = self.activ_dof[idof]
#   Check whether this dof is under Dirichlet, P.B.C., or free of B.C.
#            if (mapped_dof == self.indx_inactiv): # Dirichlet B.C.
#                # assemble global force vector contributed by Dirichlet B.C.
#                for jindx,jdof in enumerate(el_dof):
#                    jmapped_dof = self.activ_dof[jdof]
#
#                    if (jmapped_dof != self.indx_inactiv):
#                        self.Fg[self.glob_dof[jmapped_dof]] += (k_el[jindx,indx]*self.disp[Ityp(idof/self.sdim),idof%self.sdim])
            if (mapped_dof == idof): # free to assemble K. Free of constraints
                # Contribution from f_int
                self.Fg[self.glob_dof[mapped_dof]] += f_el[indx] 
                for jindx,jdof in enumerate(el_dof):
                    jmapped_dof = self.activ_dof[jdof]                    
                    if (jmapped_dof != self.indx_inactiv):
                        self.Kg[self.glob_dof[mapped_dof],self.glob_dof[jmapped_dof]] += k_el[indx,jindx]
            elif (mapped_dof != self.indx_inactiv): # rest of all the cases are periodic boundary condition
                # Contribution from f_int
                self.Fg[self.glob_dof[mapped_dof]] += f_el[indx]
                # assemble Kg and Fg
                for jindx,jdof in enumerate(el_dof):
                    jmapped_dof = self.activ_dof[jdof]
                    if (jmapped_dof != self.indx_inactiv):
                        self.Kg[self.glob_dof[mapped_dof],self.glob_dof[jmapped_dof]] += k_el[indx,jindx]
                        self.Fg[self.glob_dof[jmapped_dof]] += (k_el[jindx,indx]*self.bc.bcs['p'][str(Ityp(idof/self.sdim))][idof%self.sdim])        


#   calculate k_geom. Geometric stiffness tensor
#   Sauer's code
#        % Geometrical Stiffness
#        dNm1 = kron(dN(:,1)*dN(:,1)',I2) ;
#        dNm2 = kron(dN(:,2)*dN(:,2)',I2) ;
#        dNm3 = kron(dN(:,1)*dN(:,2)',I2) ;
#        kgeo = kgeo + ( dNm1*sig(1)+dNm2*sig(2)+(dNm3+dNm3')*sig(4) ) * detJ * wg ;
                        
    def calc_k_geom(self,dNdXs,S,w,detJ,k_geom):
        sdim = dNdXs.shape[0]
        integrator = w*detJ
        for indx, dNi in enumerate(dNdXs):
            for jindx, dNj in enumerate(dNdXs):
                k_geom[:,:] += integrator*S[indx,jindx]*np.kron(np.tensordot(dNi,dNj,axes=0),np.identity(sdim))
                    
#   loop over hex8 elements
    def loop_hex8(self):
        def assemble_B_L(F,dNdXs,B_L):
            isqrt2 = 1.0/math.sqrt(2.0)
            for i in range(8):
                B_L[0,i*3:i*3+3] = np.array([F[0,0]*dNdXs[0][i],F[1,0]*dNdXs[0][i],F[2,0]*dNdXs[0][i]])
                B_L[1,i*3:i*3+3] = np.array([F[0,1]*dNdXs[1][i],F[1,1]*dNdXs[1][i],F[2,1]*dNdXs[1][i]])
                B_L[2,i*3:i*3+3] = np.array([F[0,2]*dNdXs[2][i],F[1,2]*dNdXs[2][i],F[2,2]*dNdXs[2][i]])
                B_L[3,i*3:i*3+3] = isqrt2*np.array([F[0,0]*dNdXs[1][i]+F[0,1]*dNdXs[0][i],F[1,0]*dNdXs[1][i]+F[1,1]*dNdXs[0][i],F[2,0]*dNdXs[1][i]+F[2,1]*dNdXs[0][i]])
                B_L[4,i*3:i*3+3] = isqrt2*np.array([F[0,1]*dNdXs[2][i]+F[0,2]*dNdXs[1][i],F[1,1]*dNdXs[2][i]+F[1,2]*dNdXs[1][i],F[2,1]*dNdXs[2][i]+F[2,2]*dNdXs[1][i]])
                B_L[5,i*3:i*3+3] = isqrt2*np.array([F[0,2]*dNdXs[0][i]+F[0,0]*dNdXs[2][i],F[1,2]*dNdXs[0][i]+F[1,0]*dNdXs[2][i],F[2,2]*dNdXs[0][i]+F[2,0]*dNdXs[2][i]])
            
    #   temporary arrays            
        elem_coords = np.zeros((self.GQs['hex8'].nen,self.GQs['hex8'].sdim),dtype=Dtyp)
        conn = np.zeros((self.GQs['hex8'].nen),dtype=Ityp)
        dSdE = np.zeros((6,6),dtype=Dtyp)
        k_elm = np.zeros((self.GQs['hex8'].nen*self.GQs['hex8'].sdim,self.GQs['hex8'].nen*self.GQs['hex8'].sdim),dtype=Dtyp)
        k_geom = np.zeros((self.GQs['hex8'].nen*self.GQs['hex8'].sdim,self.GQs['hex8'].nen*self.GQs['hex8'].sdim),dtype=Dtyp)
        f_e = np.zeros((self.GQs['hex8'].nen*self.GQs['hex8'].sdim,1),dtype=Dtyp)
        f_int = np.zeros((self.GQs['hex8'].nen*self.GQs['hex8'].sdim,1),dtype=Dtyp)
        B_L = np.zeros((6,self.sdim*self.GQs['hex8'].nen),dtype=Dtyp)
#   Jacobian at each gauss point            
        for indx in range(self.nhex8):                
            self.mesh.elem_node_coords(self.mesh.hex8[indx],elem_coords)
            conn[:] = self.mesh.hex8[indx]
            self.GQs['hex8'].calc_dNdXs(elem_coords)
#   add displacement to elem_coords: elem_coords x = X + u
            for i_tmp in range(elem_coords.shape[0]):
                elem_coords[i_tmp] = elem_coords[i_tmp] + self.disp[conn[i_tmp]] 
            k_elm.fill(0.0) # initialize element sitffness matrix
            k_geom.fill(0.0)
            f_e.fill(0.0)
            f_int.fill(0.0)
            #   loop over gauss points to calculate F, S 
            for igq in range(self.GQs['hex8'].ngq):
                dSdE.fill(0.0)
                # calculate (F=dx/dX)^T 
                np.matmul(self.GQs['hex8'].dNdXs[igq],elem_coords,self.F['hex8'][indx][igq])
                self.F['hex8'][indx][igq] = self.F['hex8'][indx][igq].transpose()
                # call appropriate material models
                mat_indx = self.el_params['hex8'][indx][igq][0]
                if (mat_indx == 0): # Neo-Hookean
                    MBE_mats.Neo_Hook(dSdE,self.F['hex8'][indx][igq],self.S['hex8'][indx][igq],self.el_params['hex8'][indx][igq][1:],finite_diff=True) # only pass params
#                    print (dSdE)
                elif (mat_indx == 1): # Linear-elasticity
                    MBE_mats.lin_elasticity(dSdE,self.F['hex8'][indx][igq],self.S['hex8'][indx][igq],self.el_params['hex8'][indx][igq][1:])
                elif (mat_indx == 2): # Hook
                    MBE_mats.hook(dSdE,self.F['hex8'][indx][igq],self.S['hex8'][indx][igq],self.el_params['hex8'][indx][igq][1:])
                    
                # Assemble material stiffness matrix
                assemble_B_L(self.F['hex8'][indx][igq],self.GQs['hex8'].dNdXs[igq],B_L)
#                    print (B_L)
                self.calc_k_elm(B_L,dSdE,self.GQs['hex8'].w[igq],self.GQs['hex8'].detJ[igq],k_elm)
                self.calc_fint(B_L,self.S['hex8'][indx][igq],self.GQs['hex8'].w[igq],self.GQs['hex8'].detJ[igq],f_int)
                if (self.nlgeom):
                    self.calc_k_geom(self.GQs['hex8'].dNdXs[igq],self.S['hex8'][indx][igq],self.GQs['hex8'].w[igq],self.GQs['hex8'].detJ[igq],k_geom)                    
            if (self.nlgeom): k_elm[:,:] += k_geom[:,:]
            f_e[:] += f_int[:]
            self.assemble_Kg_Fg(conn,k_elm,f_e)

#   solve Kg*solution = Fg
#   calculate solution = inv(Kg)*Fg                    
    def solve(self):
        if (self.mat_type=="dense"):
            self.solution = np.linalg.solve(self.Kg,-self.Fg)

#   update displacement from solution vector    
    def update_disp(self):            
        for indx,dof in enumerate(self.inv_glob_dof):
            node = Ityp(dof/self.sdim)
            direction = Ityp(dof % self.sdim)
            self.disp[node,direction] += self.solution[indx]
#   Now, take nodes on periodic B.C.s into account
        for key,value in self.bc.bcs['p'].items():
            n_p = Ityp(key) # node under P.B.C.
            n_c = value[0] # corresponding node            
            for i in range(self.sdim):
                self.disp[n_p,i] += self.solution[self.glob_dof[n_c*self.sdim + i]]  
                    
    def NR_iteration(self,**kwargs):
        tol = 1.0e-7
        max_it = 20
        if (kwargs is not None):
            for key, value in kwargs.items():
                if (key == "tol"):
                    max_it = tol
                if (key == "max_it"):
                    max_it = value
        for i in range(max_it):   
            self.loop_elem()
            self.solve()
            self.update_disp()    
            resi = np.linalg.norm(self.solution) / len(self.solution)
            print("Step: "+str(i))
            print("Average norm of the increment: "+str(resi))
            if (resi< tol):
                print("Convergence attained")
                break
            else:
                continue
        if (resi>tol):
            print("Maximum iteration "+str(max_it)+" is reached, but NR did not converge.")
            print("residual = "+str(resi))
            print("tolerance = "+str(tol))


            
                
            

    def __del__(self):
        class_name = self.__class__.__name__
        print (class_name, "destroyed")
        




