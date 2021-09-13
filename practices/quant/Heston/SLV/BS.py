from scipy import optimize
import numpy as np


def CND(X):
 
    (a1,a2,a3,a4,a5) = (0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
    
    L = abs(X)
    K = 1.0 / (1.0 + 0.2316419 * L)
    w = 1.0 - 1.0 / np.sqrt(2*np.pi)*np.exp(-L*L/2.) * (a1*K + a2*K*K + a3*pow(K,3) + a4*pow(K,4) + a5*pow(K,5))
    
    if X<0:
        w = 1.0-w
 
    return w

def BlackScholes(v,CallPutFlag = 'c',S = 100.,X = 100.,T = 1.,r = 0.01):
 
    d1 = (np.log(S/X)+(r+v*v/2.)*T)/(v*np.sqrt(T))
    d2 = d1-v*np.sqrt(T)
    
    if CallPutFlag=='c':
        return S*CND(d1)-X*np.exp(-r*T)*CND(d2)
 
    else:
        return X*np.exp(-r*T)*CND(-d2)-S*CND(-d1)
    
def calc_impl_vol(price = 5., right = 'c', underlying = 100., strike = 100., time = 1., rf = 0.0):
    
    def f(x):
        out = (BlackScholes(x,CallPutFlag=right,S=underlying,X=strike,T=time,r=rf)-price)**2 
        if x < 0.0:
            out += 1000.0 * (x)**2
        return out
  
    return optimize.minimize(f,x0=0.5, tol=1e-8, method='Nelder-Mead')