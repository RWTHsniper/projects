import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import pandas as pd


def Heston_feller_cond(kappa, theta, sigma):
    return 2*kappa*theta > sigma**2

"""
Compute E[e^{1j*u*x}]
x: log(S)
u: variable for the chf
v0: v0
r: risk-free rate
rho: correlation
ka: kappa
th: theta
sig: sigma
"""
def get_Heston_chf(x, u, v0, r, rho, ka, th, sig, T):
    gam = np.sqrt(sig**2*u*(u+1j) + (ka-1j*rho*sig*u)**2)
    half_gam_T = gam*T/2.0
    coth = 1.0/np.tanh(half_gam_T)
    log_numerator = ka*th*T*(ka-1j*rho*sig*u)/sig**2 + 1j*u*(T*r+x) - (u*(u+1j)*v0)/(gam*coth + ka - 1j*rho*sig*u)
    numerator = np.exp(log_numerator)
    denominator = (np.cosh(half_gam_T) + (ka-1j*rho*sig*u)/gam*np.sinh(half_gam_T))**(2.0*ka*th/sig**2)
    return numerator/denominator

def get_damped_Heston_chf(x, v, v0, r, rho, ka, th, sig, T, alp=0.75):
    chf_inp = v - (alp+1.0)*1j
    numerator = np.exp(-r*T) * get_Heston_chf(x, chf_inp, v0, r, rho, ka, th, sig, T)
    denominator = alp**2 + alp - v**2 +1j*(2.0*alp +1.0)*v
    return numerator / denominator

def HestonMC (spot, v0, rho, kappa, theta, sigma, T, num_steps, r=0.0, return_log = True):
    # Generate a path
    vt    = np.zeros(num_steps)
    vt[0] = v0
    logSt = np.zeros(num_steps)
    logSt[0] = np.log(spot)
    dt = T / num_steps

    # Milstein scheme for volatility
    for t in range(1,num_steps):
        # Generate random Brownian Motion
        dW_indep = np.random.normal(0.0,1.0,2)
        dW_v = dW_indep[0]
        dW_logS = rho*dW_indep[0] + np.sqrt(1-rho**2)*dW_indep[1] 
        vt[t] = vt[t-1] + kappa*(theta-vt[t-1])*dt + sigma* np.sqrt(vt[t-1]*dt)*dW_v + sigma**2/4.*dt*(dW_v**2-1.)
        if vt[t] < 0.0:
            vt[t] = 0.0
        logSt[t] = logSt[t-1] + (r - vt[t-1]/2.)*dt + np.sqrt(vt[t-1]*dt)*dW_logS

    if return_log:
        return logSt, vt
    else:        
        St= np.exp(logSt)
        return St, vt
    
def get_Heston_fft_call(k, x, v0, r, rho, kappa, theta, sigma, T, alp = 0.75, N=2**12, dk=0.01):# 2048*5
    u_list = np.array(range(N-1))
    # be careful about scales of dk and dv
#    dk = np.abs(k)*0.01
    dv = 2.0*np.pi/(N*dk)  # N*dv = 2*pi/dk

#    dk = 2.0*np.pi/(N*dv) # nu nu*zeta = 2*pi/N
    
#    print(dv*(N-1))
#    print('dk ds ', dk, ' ',dv)
    k_list = np.array([x + dk*u_elem - N*dk/2.0 for u_elem in u_list])
    v_list = np.array(range(N-1)) * dv
#    print('maxv ', v_list[-1])
#    if np.isinf(np.sinh(v_list[-1]**2)) | np.isinf(np.cosh(v_list[-1]**2)) | np.isinf(np.exp(v_list[-1]**2)):
#        print('Infinite ',v_list[-1]**2)
#        raise Exception('inifite')

    x_list = [] # values in frequency domain
    for (j,vj) in enumerate(v_list): # j=0 to N-1
        if j==0:
            kroneker_delta = 1.0
        else:
            kroneker_delta = 0.0        
        simpson_coeff = dv/3.0*(3.0 + (-1)**j - kroneker_delta)
        xj = np.exp(1j*(N/2.0*dk-x)*vj)*simpson_coeff*get_damped_Heston_chf(x, vj, v0, r, rho, kappa, theta, sigma, T, alp)
        x_list.append(xj)

    fft_prices = np.exp(-alp*k)/np.pi*fft(x_list).real
#    price_index = np.where(k_list == k)[0]
#    fft_price = fft_prices[int(N/2)]
    # interpolate results
    fft_price = np.interp(np.exp(k), np.exp(k_list), fft_prices)
#    print(fft_price)
    return fft_price