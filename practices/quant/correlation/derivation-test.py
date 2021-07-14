#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy
sp = sympy
import numpy as np
from scipy.integrate import odeint,RK45,solve_ivp
import matplotlib.pyplot as plt


# In[12]:


def get_E_vAvB(t,alp,beta,sigma,v0):
    res = 1.0
    for i in range(len(beta)):
        res *= sp.exp(-beta[i]*t)*(v0[i] - alp[i]/beta[i]) + alp[i]/beta[i]
    return res

def Riccati_solution(a,b,c,tau,init_val):
    g = sp.sqrt(b**2 - 4.0*a*c)
    ex = sp.exp(g*tau) - 1
    num = 2*g*init_val + ((b+g)*init_val+2*c)*ex
    denom = 2*g - (2*a*init_val+b-g)*ex
    return num/denom

def Riccati_solution_vec(a,b,c,tau,init_val):
    n = len(a)
    res = np.zeros(n)
    for i in range(n):
        res[i] = Riccati_solution(a[i],b[i],c[i],tau,init_val[i])
    return res
    
def int_Riccati_solution(a,b,c,tau,ricc_init_val):
    g = sp.sqrt(b**2 - 4.0*a*c)
    ex = sp.exp(g*tau) - 1
    term1 = (-b+g)*tau/(2.0*a) 
    term2 = sp.log(2*g/(2*g+(g-b-2*a*ricc_init_val)*ex))/a
#     term2 = sp.log(2*g/sp.Abs(2*g+(g-b-2*a*ricc_init_val)*ex))/a
    return term1 + term2

def int_Riccati_solution_vec(a,b,c,tau,ricc_init_val):
    n = len(a)
    res = np.zeros(n)
    for i in range(n):
#        print(int_Riccati_solution(a[i],b[i],c[i],tau,ricc_init_val[i]))
        res[i] = int_Riccati_solution(a[i],b[i],c[i],tau,ricc_init_val[i])
    return res


# In[26]:


params = {
# chf params    
    'tau': 1, # investment time
    'h': 60, # horizon in the future
    'w': [1,1,1],
    'u': [1,1],
    'psi': [1,1],
# treasury params (US)
    'l0': -0.127061, # constant shift
#     've': np.array([0.2715618,0.0195524,0.0009720]),
    'th': np.array([0.04783331, 0.07757864, 0.00661224]),
    'ka': np.array([5.6772530,0.2520333,0.147]),
    'si': np.array([0.0181427,0.0422960,0.034]), # adjusted si[2]
    'x0': np.array([0.05095958,0.06725220,0.00961570]),
# Equity parameters (US LC, OS)
    's0': np.array([0.0,0.0]), # initial value
    'mu0': np.array([0.071,0.075]), # official value
    'mu1': np.array([0.066,0.055]) + np.array([-0.5,-0.5]), # I think just -0.5 is too strict? High volatility -> Always stock crush
    # mu1 = np.array([-0.001,-0.002]) + np.array([-0.5,-0.5]) # I think just -0.5 is too strict?
    'v0': np.array([0.018,0.013]),
    'alpha': np.array([0.021,0.015]),
    'beta': np.array([1.279,1.115]),
    'sigma': np.array([0.232,0.099]), # volatility of variance
    'sAsB': 0.00728,
    'sAvA': -0.444, 
    'sAvB': -0.406,
    'sBvA': -0.345,
    'sBvB': -0.421,
    'vAvB': 0.903,
}
names = ['th','ka','si','x0','s0','mu0','mu1','v0','alpha','beta','sigma','u','w','psi']
for name in names:
    for i in range(len(params[name])):
        params[name+str(i)] = params[name][i]
    

# # Characteristic Function Derivation

# In[3]:


a, b, c, tau, f0 = sp.symbols("a b c tau f0")
R_sol = Riccati_solution(a,b,c,tau,f0)
int_R_sol = int_Riccati_solution(a,b,c,tau,f0)


# In[4]:


# u = np.array([sp.symbols('u0'),sp.symbols('u1')])
# w = np.array([sp.symbols('w0'),sp.symbols('w1'),sp.symbols('u2')])
# psi = np.array([sp.symbols('psi0'),sp.symbols('psi1')])
u0,u1 = sp.symbols("u0,u1") 
u = np.array([u0,u1])
w0,w1,w2 = sp.symbols("w0,w1,w2")
w = np.array([w0,w1,w2])
psi0,psi1 = sp.symbols("psi0,psi1")
psi = np.array([psi0,psi1])
tau, h, l0 = sp.symbols('tau, h, l0')
# treasury
ve = np.array([sp.symbols('ve0'),sp.symbols('ve1'),sp.symbols('ve2')])
ka = np.array([sp.symbols('ka0'),sp.symbols('ka1'),sp.symbols('ka2')])
th = np.array([sp.symbols('th0'),sp.symbols('th1'),sp.symbols('th2')])
si = np.array([sp.symbols('si0'),sp.symbols('si1'),sp.symbols('si2')])
x0 = np.array([sp.symbols('x00'),sp.symbols('x01'),sp.symbols('x02')])
# ve, ka, si, x0 = sp.symbols("ve, ka, si, x0")
# Presume one factor 
# equity
s0 = np.array([sp.symbols('s00'),sp.symbols('s01')])
mu0 = np.array([sp.symbols('mu00'),sp.symbols('mu01')])
mu1 = np.array([sp.symbols('mu10'),sp.symbols('mu11')])
v0 = np.array([sp.symbols('v00'),sp.symbols('v01')])
alpha = np.array([sp.symbols('alpha0'),sp.symbols('alpha1')])
beta = np.array([sp.symbols('beta0'),sp.symbols('beta1')])
sigma = np.array([sp.symbols('sigma0'),sp.symbols('sigma1')])
correlation = np.matrix([[1,sp.symbols('sAsB'),sp.symbols('sAvA'),sp.symbols('sAvB')],
                        [sp.symbols('sAsB'),1,sp.symbols('sBvA'),sp.symbols('sBvB')],
                        [sp.symbols('sAvA'),sp.symbols('sBvA'),1,sp.symbols('vAvB')],
                        [sp.symbols('sAvB'),sp.symbols('sBvB'),sp.symbols('vAvB'),1]]) # uA uB vA vB

# Check correlation
# In[13]:


B = []
for i in range(3):
    B.append(R_sol.subs({a:0.5*si[i]**2,b:-ka[i],c:u.sum(),f0:w[i]}))
B = np.array(B)
C = []
sAsB = correlation[0,1]
sAvA = correlation[0,2]
sAvB = correlation[0,3]
sBvA = correlation[1,2]
sBvB = correlation[1,3]
i=0
C.append(R_sol.subs({a:0.5*sigma[i]**2,b:u[i]*sAvA*sigma[i]-beta[i],c:0.5*u[i]**2+u[i]*mu1[i],f0:psi[i]}))
i=1
C.append(R_sol.subs({a:0.5*sigma[i]**2,b:u[i]*sBvB*sigma[i]-beta[i],c:0.5*u[i]**2+u[i]*mu1[i],f0:psi[i]}))
# E_vAvB = sp.symbols('E_vAvB')
E_vAvB = get_E_vAvB(tau,alpha,beta,sigma,v0)
# E_vAvB = (alpha/beta).prod()
# A = (((u*mu0).sum() + (l0*u).sum() + E_vAvB*u.prod()*sAsB)*tau).sum() # check
A = ((u*mu0).sum() + (l0*u).sum() + E_vAvB*u.prod()*sAsB)*tau
int_B = []
int_C = []
i=0
int_B.append(int_R_sol.subs({a:0.5*si[i]**2,b:-ka[i],c:u.sum(),f0:w[i]}))
int_C.append(int_R_sol.subs({a:0.5*sigma[i]**2,b:u[i]*sAvA*sigma[i]-beta[i],c:0.5*u[i]**2+u[i]*mu1[i],f0:psi[i]}))
i=1
int_B.append(int_R_sol.subs({a:0.5*si[i]**2,b:-ka[i],c:u.sum(),f0:w[i]}))
int_C.append(int_R_sol.subs({a:0.5*sigma[i]**2,b:u[i]*sBvB*sigma[i]-beta[i],c:0.5*u[i]**2+u[i]*mu1[i],f0:psi[i]}))
i=2
int_B.append(int_R_sol.subs({a:0.5*si[i]**2,b:-ka[i],c:u.sum(),f0:w[i]}))
for ind,elem in enumerate(B):
    A += ka[ind]*th[ind]*int_B[ind]
i=0
A += (alpha[i] + u[1]*sBvA*sigma[i]*E_vAvB)*int_C[i] # check
i=1
A += (alpha[i] + u[0]*sAvB*sigma[i]*E_vAvB)*int_C[i] # check
chf = sp.exp(A+np.dot(u,s0)+np.dot(B,x0)+np.dot(C,v0))
# Check implementation of A

params2 = params.copy()
params2['tau'] = 1
print(chf.subs(params2))
# # Chf for Log-Return

# In[14]:

p_logret = params.copy()
del p_logret['u0']
del p_logret['u1']
del p_logret['w0']
del p_logret['w1']
del p_logret['w2']
del p_logret['psi0']
del p_logret['psi1']
p_logret['tau'] = tau
#del p_logret['tau']


#inp = {tau:h,u0:0,u1:0,w0:B[0],w1:B[1],w2:B[2],psi0:C[0],psi1:C[1]}

## First affine form

chf_params ={w0:0,w1:0,w2:0,psi0:0,psi1:0}
A0 = A.subs(chf_params)
B0 = []
for elem in B:
    B0.append(elem.subs(chf_params))
C0 = []
for elem in C:
    C0.append(elem.subs(chf_params))
B0 = np.array(B0)
C0 = np.array(C0)


## Second affine form
chf_params2 = {tau:h,
               u0:0,u1:0,
               w0:B0[0],w1:B0[1],w2:B0[2],
               psi0:C0[0],psi1:C0[1]}

A2 = A.subs(chf_params2)
B2 = []
for elem in B:
    B2.append(elem.subs(chf_params2))
C2 = []
for elem in C:
    C2.append(elem.subs(chf_params2))
B2 = np.array(B2)
C2 = np.array(C2)


chf_logret = sp.exp(A0+A2+np.dot(B2,x0)+np.dot(C2,v0))

zero_args = {u0:0,u1:0,w0:0,w1:0,w2:0,psi0:0,psi1:0,tau:1}
p_logret['h'] = h
sim_chf_logret = chf_logret.subs(p_logret) # insert model parameters
sim_chf_u0 = sim_chf_logret.diff(u0)
sim_chf_u1 = sim_chf_logret.diff(u1)
sim_chf_u0u1 = sim_chf_u0.diff(u1)
mean_ret_u0 = sim_chf_u0.subs(zero_args)
mean_ret_u1 = sim_chf_u1.subs(zero_args)
mean_ret_u0u1 = sim_chf_u0u1.subs(zero_args)


time = np.linspace(0.0,30.0,num=31)
cov = []

for hor in time:
    inp = {h:hor}
    ret_cov = mean_ret_u0u1.subs(inp)-mean_ret_u0.subs(inp)*mean_ret_u1.subs(inp)
    cov.append(ret_cov)

plt.plot(cov)
plt.title('Theoretical Covariance')

# In[16]:

# validation of sim_chf_logret

time = np.linspace(0.0,30.0,num=31)
ret = []
for hor in time:
    chf_ret =  sim_chf_logret.subs({h:hor,u0:1,u1:1,w0:1,w1:1,w2:1,psi0:1,psi1:1,tau:1})
    ret.append(chf_ret)

plt.figure()
plt.plot(ret)
plt.title('Theoretical chf for return')



    
