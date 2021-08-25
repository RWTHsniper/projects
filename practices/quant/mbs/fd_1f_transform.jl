
#include("/Users/jungjaeyong/projects/practices/quant/mbs/fd_1f.jl")

import SparseArrays; sp = SparseArrays
import IterativeSolvers; is = IterativeSolvers
import LinearAlgebra; la = LinearAlgebra
import Interpolations; ip = Interpolations


# CIR process functions

function r_a_aux0(dl,th,ka,si,tau,v=1.0,u=0.0)
    g_aux = sqrt(ka^2 + 2*dl*v*si^2)
    ex = exp(g_aux*tau) - 1.0
    c = 2.0*g_aux + ex*(ka + g_aux - u*si^2)
    return ka*th*(tau*(ka + g_aux) + 2*log(2*g_aux/c))/(si^2)
end

function r_a_aux0(dl::Vector,th::Vector,ka::Vector,si::Vector,tau,v=1.0,u=0.0)
    res = zeros(size(dl))
    for i in 1:length(res)
        res[i] = r_a_aux0(dl[i],th[i],ka[i],si[i],tau,v,u)
    end
    return res
end

function r_b_aux0(dl,ka,si,tau,v=1.0,u=0.0)
    g_aux = sqrt(ka^2 + 2*dl*v*si^2)
    ex = exp(g_aux*tau) - 1.0
    c = 2.0*g_aux + ex*(ka + g_aux - u*si^2)
    return (2*u*g_aux - ex*(2*dl*v + u*(ka - g_aux)))/c
end

function r_b_aux0(dl::Vector,ka::Vector,si::Vector,tau,v=1.0,u=0.0)
    res = zeros(size(dl))
    for i in 1:length(res)
        res[i] = r_b_aux0(dl[i],ka[i],si[i],tau,v,u)
    end
    return res
end

function get_n_nonzeros(Nx)
    if typeof(Nx) <: Real 
        res = (Nx-2)*3 + 2*2
    else
        res = 0
        dim = length(Nx)
        if dim == 3
            res += prod(Nx .- 2)*7 # inner nodes
            res += 2^dim*4 # vertex nodes
            res += 2*(dim-1)*(sum(Nx) - 2*dim)*5 # edges
            res += 2*((Nx[1]-2)*(Nx[2]-2)+(Nx[1]-2)*(Nx[3]-2)+(Nx[2]-2)*(Nx[3]-2))*6 # outer surface
        end
    end
    return res
end

function get_x(x_min,dx,ind)
    return x_min + dx*(ind-1)
end



"""
Ic: I for the matrix for the current step
Jc: J for the matrix for the current step
Vc: Val for the matrix for the current step
"""
function get_trans_mat_COO!(Nx, params, dx,dt,I,J,Val,Ic,Jc,Vc,x_min)
    si = params[:si]
    ve = params[:ve]
    ka = params[:ka]
    k = params[:k] # strike interest rate
    gamma = params[:gamma]
    count = 1 # counter to point an element in I,J,Val
    # first index
    ind = 1 
    indm = ind - 1
    indp = ind + 1  
    x_ind = get_x(x_min,dx,ind)
    if k > x_ind
        C2 = (ve - (gamma + params[:hb]*si^2)*x_ind)/dx # Qx: (ve - ka*x)/dx, ve/dx - ka*(x_min)/dx
    else
        C2 = (ve - (gamma + params[:hb1]*si^2)*x_ind)/dx # Qx: (ve - ka*x)/dx
    end
    # C2 = ve/dx - ka*(x_min)/dx # Qx: (ve - ka*x)/dx
    # C3 = (1.0-gamma)*(x_min) + gamma*k # Q
    C3 = 0.0
    I[count] = ind; J[count] = ind;
    Val[count] = 2.0/dt + C2 + C3
    center = count
    count += 1
    I[count] = ind; J[count] = indp; Val[count] = -C2
    count += 1
    for ind in 2:(Nx-1)
        x_ind = get_x(x_min,dx,ind)
        if k > x_ind  
            C2 = (ve - (gamma +params[:hb]*si^2)*x_ind)/(2.0*dx) # ve/(2.0*dx) - ka*((ind-1)*dx+x_min)/(2.0*dx)
        else
            C2 = (ve - (gamma + params[:hb1]*si^2)*x_ind)/(2.0*dx) # ve/(2.0*dx) - ka*((ind-1)*dx+x_min)/(2.0*dx)
        end
        C1 = si^2 *(x_ind)/(2.0*dx^2)
        C3 = 0.0
        indm = ind - 1
        indp = ind + 1
        I[count] = ind; J[count] = indm
        Val[count] = -C1 + C2
        count += 1
        I[count] = ind; J[count] = ind
        Val[count] = 2.0*C1 + C3 + 2.0/dt
        count += 1
        I[count] = ind; J[count] = indp
        Val[count] = -C1 - C2
        count += 1
    end
    ind = Nx
    indm = ind - 1
    indp = ind + 1
    x_ind = get_x(x_min,dx,ind)
    if k > x_ind
        C2 = (ve - (gamma + params[:hb]*si^2)*x_ind)/dx # Qx: (ve - ka*x)/dx, ve/dx - ka*(x_min)/dx
    else
        C2 = (ve - (gamma + params[:hb1]*si^2)*x_ind)/dx # Qx: (ve - ka*x)/dx
    end    
    C3 = 0.0
    I[count] = ind; J[count] = indm
    Val[count] = C2
    count += 1
    I[count] = ind; J[count] = ind
    Val[count] = 2.0/dt - C2 + C3
    center = count
    count += 1

    # matrix for current step
    Ic .= I
    Jc .= J
    Vc .= -Val
    # corrections on centeral values
    for i in 1:length(I)
        if I[i] == J[i] # node center
            Vc[i] += 4.0/dt # 2.0*2.0
        end
    end

end



"""
Ic: I for the matrix for the current step
Jc: J for the matrix for the current step
Vc: Val for the matrix for the current step
"""
function get_mat_COO!(Nx, params, dx,dt,I,J,Val,Ic,Jc,Vc,x_min)
    si = params[:si]
    ve = params[:ve]
    ka = params[:ka]
    k = params[:k] # strike interest rate
    count = 1 # counter to point an element in I,J,Val
    # first index
    ind = 1 
    indm = ind - 1
    indp = ind + 1  
    x_ind = x_min
    if k > x_ind
        gamma = params[:gamma]
    else
        gamma = 0.0
    end
    C2 = ve/dx - ka*(x_min)/dx
    C3 = (1.0-gamma)*(x_min) + gamma*k
    I[count] = ind; J[count] = ind;
    Val[count] = 2.0/dt + C2 + C3
    center = count
    count += 1
    I[count] = ind; J[count] = indp; Val[count] = -C2
    count += 1
    for ind in 2:(Nx-1)
        x_ind = x_min + (ind-1)*dx
        if k > x_ind  
            gamma = params[:gamma]
        else
            gamma = 0.0
        end
        C1 = si^2 *(x_ind)/(2.0*dx^2)
        C2 = ve/(2.0*dx) - ka*((ind-1)*dx+x_min)/(2.0*dx)
        C3 = (1.0-gamma)*(x_ind) + gamma*k
        indm = ind - 1
        indp = ind + 1
        I[count] = ind; J[count] = indm
        Val[count] = -C1 + C2
        count += 1
        I[count] = ind; J[count] = ind
        Val[count] = 2.0*C1 + C3 + 2.0/dt
        count += 1
        I[count] = ind; J[count] = indp
        Val[count] = -C1 - C2
        count += 1
    end
    ind = Nx
    indm = ind - 1
    indp = ind + 1
    x_ind = x_min + (ind-1)*dx
    if k > x_ind
        gamma = params[:gamma]
    else
        gamma = 0.0
    end
    C2 = ve/dx - ka*(x_ind)/dx
    C3 = (1.0-gamma)*(x_ind) + gamma*k
    I[count] = ind; J[count] = indm
    Val[count] = C2
    count += 1
    I[count] = ind; J[count] = ind
    Val[count] = 2.0/dt - C2 + C3
    center = count
    count += 1

    # matrix for current step
    Ic .= I
    Jc .= J
    Vc .= -Val
    # corrections on centeral values
    for i in 1:length(I)
        if I[i] == J[i] # node center
            Vc[i] += 4.0/dt # 2.0*2.0
        end
    end

end


function compute_fd!(Vc,Vn,matc,mat,t_space)
    tmp = zeros(size(Vc))
    for (ind, tau) in enumerate(t_space)
        tmp .= matc*Vc
        # Vn .= mat\tmp
        is.gmres!(Vn,mat,tmp)
        Vc .= Vn
    end
    return Vn
end

function get_quadratic_sol(a,b,c)
    res = ((-b + sqrt(b^2 -4*a*c+0im))/(2.0*a),(-b - sqrt(b^2 -4*a*c+0im))/(2.0*a))
    return res
end


#=
Think about how to contruct an FD scheme
at tau=0, initial conditions are given.
At low x and upper x, I impose zero-gamma BC


_______________________________ x = x_max
|                             |
|                             |
|                             |
|                             |
_______________________________ x = x_min
tau = 0                       tau = T

In Julia, it is column-major indexing. 
Since I am using a COO sparse matrix, it does not really matter. 
Otherwise, I would carefully choose CSR or CSC.
Therefore, let's start with COO first.
=#

params = Dict(
# contract params
    :T => 30, # test
# Treasury
    :l0 => 0.0,
    :ve => 0.0027646771594520936,
    :ka => 0.06491552810007564, 
    :si => 0.03362592979733442, 
    :x0 => 0.000100914, # test
    :x0 => 0.02, # test
# harzard rate process (PSA params)
# ho(t) parameters
    :a => 0.0,
    :b => 0.0, 
    :gamma => 20.0,
    :k => 0.02, # prepayment strike
    # :k => 0.0, # prepayment strike
    :T_asterisk => 1000.0, # prepayment date test # test
    )

th = params[:ve]/params[:ka]
k = params[:k]
x_min = 0.0000001 # CIR
x_max = 0.6
x0 = params[:x0]
x_min = x0/10.0 # CIR
x_min = 1e-6
# x_max = k*10.0
grid_d = Dict(
    :Nt => 30*12, # num of timesteps
    # :Nx => [3,3,3], # num of state variables
    :Nx => 64, # num of state variables
    :Nx => 256, # num of state variables
    # :Nx => 256, # num of state variables
    # :Nx => 1024, # num of state variables
    # :Nx => [128,128,32], # num of state variables
    :x_min => x_min,
    :x_max => x_max,
)

H_params = [0.5*params[:si]^2,params[:gamma],-1.0+params[:gamma]]
hb = get_quadratic_sol(H_params...)
hb_ind = argmin(abs.(hb))
# hb_ind = argmax(abs.(hb))
hb = hb[hb_ind] # b in H transform
H_params_1 = [0.5*params[:si]^2,params[:gamma],-1.0]
hb1 = get_quadratic_sol(H_params_1...)
hb1_ind = argmin(abs.(hb1))
hb1 = hb1[hb1_ind] # b1 in H1 transform
t_params = copy(params)
t_params[:hb] = hb
t_params[:hb1] = hb1

dx = (grid_d[:x_max] - grid_d[:x_min]) / (grid_d[:Nx] - 1.0)
dt = params[:T] / grid_d[:Nt]
Nx = grid_d[:Nx]
num_nonzeros = get_n_nonzeros(Nx)
I = zeros(num_nonzeros)
J = zeros(num_nonzeros)
Val = zeros(num_nonzeros)
Ic = zeros(num_nonzeros)
Jc = zeros(num_nonzeros)
Vc = zeros(num_nonzeros)
# get_mat_COO!(Nx, params, dx,dt,I,J,Val,Ic,Jc,Vc,x_min)
get_trans_mat_COO!(Nx, t_params, dx,dt,I,J,Val,Ic,Jc,Vc,x_min)

mat = sp.sparse(I,J,Val) # next step
matc = sp.sparse(Ic,Jc,Vc) # current step

T = params[:T]
tau_space = collect(0.0:dt:T)

# Try to solve it using an iterative IterativeSolvers

# 1. fQ(ru) = 1
fQ(ru) = 1.0
Q0 = zeros(Nx) # initial value
Qn = zeros(Nx) # next
for (ind,elem) in enumerate(Q0)
    Q0[ind] = fQ(elem)
end
# Transform to H
hQ0 = copy(Q0)
tau = 0.0
for (ind,elem) in enumerate(Q0)
    x_ind = get_x(x_min,dx,ind)
    if x_ind < k
        hQ0[ind] = Q0[ind]*exp(tau*(t_params[:gamma]*k + t_params[:hb]*t_params[:ve]) + t_params[:hb]*x_ind)
    else
        hQ0[ind] = Q0[ind]*exp(tau*(t_params[:gamma]*k + t_params[:hb1]*t_params[:ve]) + t_params[:hb1]*x_ind)
    end
end


hQc = deepcopy(hQ0) # current
hQn = zeros(Nx) # current
compute_fd!(hQc,hQn,matc,mat,@view tau_space[2:end])

# inverse transform from H to Q
T = t_params[:T]
for (ind,elem) in enumerate(hQn)
    x_ind = get_x(x_min,dx,ind)
    if x_ind < k
        Qn[ind] = hQn[ind]/exp(T*(t_params[:gamma]*k + t_params[:hb]*t_params[:ve]) + t_params[:hb]*x_ind)
    else
        Qn[ind] = hQn[ind]/exp(T*(t_params[:gamma]*k + t_params[:hb1]*t_params[:ve]) + t_params[:hb1]*x_ind)
    end
end


x0_ind = round(Int,(x0 -x_min)/dx + 1)
println("Qn at t=0 ", Qn[x0_ind])


# Analytic Q
dl = 1
ka = params[:ka]
si = params[:si]
tau = T
A = r_a_aux0(dl,th,ka,si,tau)
B =r_b_aux0(dl,ka,si,tau)
bond_price = exp(sum(A+B*x0))
println("Analytic solution for bond price: ",bond_price)

# compute R
fR(ru) = ru
R0 = zeros(Nx) # initial value
Rn = zeros(Nx) # next
for (ind,elem) in enumerate(R0)
    x_ind = x_min + dx*(ind-1)
    R0[ind] = fR(x_ind)
end
Rc = deepcopy(R0) # current
compute_fd!(Rc,Rn,matc,mat,@view tau_space[2:end])
x0_ind = round(Int,(x0 -x_min)/dx + 1)
println("Rn at t=0 ", Rn[x0_ind])
