import numpy as np
import scipy
import copy


def get_ind(Nx,pos,order='C'):
    """
    Compute the index in the global scope
    Nx: grid
    pos: position in the grid
    order: 'C': C-lang, 'F': Fortran
    """
    res = 0
    multiplier = 1
    if order == 'C':
        for ind in reversed(range(len(Nx))):
            if (pos[ind] >= Nx[ind]) | (pos[ind] < 0):
                raise("Index range is wrong")
            res += pos[ind]*multiplier
            multiplier *= Nx[ind]
    elif order == 'F':
        for ind in range(len(Nx)):
            if (pos[ind] > Nx[ind]) | (pos[ind] <= 0):
                raise("Index range is wrong")
            res += pos[ind]*multiplier
            multiplier *= Nx[ind]
    return int(res)

def get_pos(Nx,ind,order='C'):
    if order == 'C':
        if ind >= Nx.prod():
            print("index is larger than maximum")
            raise(ValueError)
        len_Nx = len(Nx)
        mod = ind
        pos = np.zeros(len_Nx)
        for i in range(len_Nx):
            divider = Nx[i+1:].prod()
            print(mod,divider)
            q, mod = divmod(mod, divider)
            print(i,q)
            pos[i] = q
    elif order == 'F':
        if ind > Nx.prod():
            print("index is larger than maximum")
            raise(ValueError)
        len_Nx = len(Nx)
        mod = ind-1
        pos = np.zeros(len_Nx)
        for i in reversed(range(len_Nx)):
            divider = Nx[0:i].prod()
            print(mod,divider)
            q, mod = divmod(mod, divider)
            print(i,q)
            pos[i] = q
        pos += 1.0
        
    return pos


def get_x(dx,pos,x_min):
    return dx*pos + x_min

def get_v(dx,pos):
    return dx[1]*pos[1]

def get_x_ind(x,dx,x_min):
    ind = round((x-x_min)/dx)
    return int(ind)

# My sparse matrix format

class tridiagonal:
    def __init__(self, diag_1=None, diag_2=None, diag_3=None, matrix=None, n=None):
        if (type(diag_1) == np.ndarray) and (type(diag_2) == np.ndarray) and (type(diag_3) == np.ndarray):
            self.n = len(diag_2)
            self.diag_1 = diag_1
            self.diag_2 = diag_2
            self.diag_3 = diag_3
        elif n is not None:
            self.n = n
            self.diag_1 = np.zeros(n-1)
            self.diag_2 = np.zeros(n)
            self.diag_3 = np.zeros(n-1)
        elif matrix is not None:
            self.n = matrix.shape[0]
            self.diag_1 = np.zeros(self.n-1) 
            self.diag_2 = np.zeros(self.n)
            self.diag_3 = np.zeros(self.n-1)
            for i in range(self.n):
                if i==0:
                    self.diag_2[i] = matrix[i,i]
                    self.diag_3[i] = matrix[i,i+1]
                elif i==self.n-1:
                    self.diag_1[i-1] = matrix[i,i-1]
                    self.diag_2[i] = matrix[i,i]        
                else:
                    self.diag_1[i-1] = matrix[i,i-1]
                    self.diag_2[i] = matrix[i,i]        
                    self.diag_3[i] = matrix[i,i+1]
    
    def __mul__(self, other):
        if type(other) == np.ndarray:
            if len(other) != self.n:
                print("Wrong length for the multiplication with a tridiagonal matrix")
            else:
                res = np.zeros(self.n)
                for i in range(self.n):
                    if i==0:
                        res[i] = self.diag_2[i]*other[i] + self.diag_3[i]*other[i+1]
                    elif i== self.n-1:
                        res[i] = self.diag_1[i-1]*other[i-1] + self.diag_2[i]*other[i]                    
                    else:
                        res[i] = self.diag_1[i-1]*other[i-1] + self.diag_2[i]*other[i] + self.diag_3[i]*other[i+1] 
        elif np.isscalar(other):
            res = tridiagonal(self.diag_1,self.diag_2,self.diag_3)
            res.diag_1 *= other
            res.diag_2 *= other
            res.diag_3 *= other            
        return res
                    
    def solve(self, rhs):
        res=scipy.linalg.lapack.dgtsv(self.diag_1,self.diag_2,self.diag_3,rhs)
        info = res[4]
        if info != 0:
            print("dgtsv solution won't be complete: ",info)
        return res[3]

    def __getitem__(self,key):
        i_row = key[0]
        i_col = key[1]
        if (i_col < i_row-1) or (i_row+1< i_col) or (len(key) != 2):
            raise(IndexError)
        if i_col < i_row:
            return self.diag_1[i_row-1]
        elif i_col == i_row:
            return self.diag_2[i_row]            
        else:
            return self.diag_3[i_row]
            
    def __setitem__(self,key,value):
        i_row = key[0]
        i_col = key[1]
        if (i_col < i_row-1) or (i_row+1< i_col) or (len(key) != 2):
            raise(IndexError)
        if i_col < i_row:
            self.diag_1[i_row-1] = value
        elif i_col == i_row:
            self.diag_2[i_row] = value            
        else:
            self.diag_3[i_row] = value

    def __add__(self, other):
        if type(other) == tridiagonal:
#             res = tridiagonal(self.diag_1,self.diag_2,self.diag_3)
            res = copy.deepcopy(self)
            res.diag_1 += other.diag_1
            res.diag_2 += other.diag_2
            res.diag_3 += other.diag_3
        return res

    def __sub__(self, other):
        if type(other) == tridiagonal:
#             res = tridiagonal(self.diag_1,self.diag_2,self.diag_3)
            res = copy.deepcopy(self)
            res.diag_1 -= other.diag_1
            res.diag_2 -= other.diag_2
            res.diag_3 -= other.diag_3
        return res
        
    def __rmul__(self, other):
        if np.isscalar(other):
#             res = tridiagonal(self.diag_1,self.diag_2,self.diag_3)
            res = copy.deepcopy(self)
            res.diag_1 *= other
            res.diag_2 *= other
            res.diag_3 *= other
        return res
            