
# include("/Users/jungjaeyong/projects/practices/quant/mbs/fd.jl")


import SparseArrays; sp = SparseArrays
import IterativeSolvers; is = IterativeSolvers


function get_ind(N, pos)
    res = 0
    # for i in 1:(length(pos)-1)
        # res += (pos[i]-1)*prod(N[i+1:end])
    # end
    # res += pos[end]
    res = pos[1]
    for i in 2:(length(pos))
        res += (pos[i]-1)*prod(N[1:i-1])
    end
    return res
end

function get_pos(Nx, ind)
    effective_dim = 0
    pos = ones(size(Nx))
    if ind <= Nx[1] 
        effective_dim = 1
        pos[1] = ind
        return pos
    elseif ind <= Nx[2]*Nx[1]
        effective_dim = 2
    else
        effective_dim = 3
    end
    # inverse function for get_ind
    for i in effective_dim:-1:2
        factor = prod(Nx[1:i-1])
        q, r = divrem(ind, factor)
        if r == 0
            pos[i] = min(Nx[i], q)
        else
            pos[i]  = min(Nx[i], q + 1)
        end
        # @show i, pos[i], factor
        ind -= factor * (pos[i]-1)
    end
    pos[1] = ind
    return pos
end

#=
# Test codes
# test get_pos and get_ind
n_dof = prod(Nx)
for i in 1:n_dof
    pos = get_pos(Nx, i)
    ind = get_ind(Nx,pos); ind = convert(Int, ind)
    if i != ind
        println("Different! ", i, " ", ind)
    end
end

# test get_n_nonzeros and get_n_nonzeros_iter
@show get_n_nonzeros(Nx) == get_n_nonzeros_iter(Nx)
=#


function compute_fd!(V0,Vc,Vn,matc,mat,t_space)
    Vc .= V0
    for (ind, tau) in enumerate(t_space)
        Vc .= matc*Vc
        Vn .= mat\Vc
        Vc .= Vn
    end
    
end

# Let's go


# Define initial condition (tau=0, t=T)
function get_n_nonzeros_iter(Nx)
    res = 0
    for c in 1:Nx[3]
        for b in 1:Nx[2]
            for a in 1:Nx[1]
                res += 1
                if a == 1 || a == Nx[1]
                    res += 1
                else
                    res += 2
                end
                if b == 1 || b == Nx[2]
                    res += 1
                else
                    res += 2
                end
                if c == 1 || c == Nx[3]
                    res += 1
                else
                    res += 2
                end
            end
        end
    end    
    return res
end

function get_n_nonzeros(Nx)
    res = 0
    dim = length(Nx)
    if dim == 3
        res += prod(Nx .- 2)*7 # inner nodes
        res += 2^dim*4 # vertex nodes
        res += 2*(dim-1)*(sum(Nx) - 2*dim)*5 # edges
        res += 2*((Nx[1]-2)*(Nx[2]-2)+(Nx[1]-2)*(Nx[3]-2)+(Nx[2]-2)*(Nx[3]-2))*6 # outer surface
    end
    return res
end



"""
Ic: I for the matrix for the current step
Jc: J for the matrix for the current step
Vc: Val for the matrix for the current step
"""
function get_mat_COO!(Nx, params, dx,dt,I,J,Val,Ic,Jc,Vc)
    si = params[:si]
    ve = params[:ve]
    ka = params[:ka]
    gamma = params[:gamma]
    k = params[:k] # strike interest rate
    count = 1 # counter to point an element in I,J,Val
    for c in 1:Nx[3]
        for b in 1:Nx[2]
            for a in 1:Nx[1]
                # try
                pos = (a,b,c)
                ind = get_ind(Nx, pos)
                # Compute vectors of C1, C2, C3
                C1 = si.^2 .*(pos.-1)./(2.0*dx)
                C2 = ve./(2.0*dx) - ka.*(pos.-1)/2.0
                C3 = (1.0-gamma).*(pos.-1).*dx .+ gamma*k
                C3[1:2] .= 0.0 # Use only the third factor for the harzard process
                # a
                p = a
                dof = 1
                indm = pos .- (1,0,0)
                indm = get_ind(Nx,indm)
                indp = pos .+ (1,0,0)
                indp = get_ind(Nx,indp)
                if 1 < p < Nx[dof] # inner
                    I[count] = ind; J[count] = indm; Val[count] = -C1[dof] + C2[dof] 
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    I[count] = ind; J[count] = ind; Val[count] = 2.0*C1[dof] + C3[dof] + 2.0/dt
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count] + 4.0/dt # offset
                    center = count # center node
                    count += 1
                    I[count] = ind; J[count] = indp; Val[count] = -C1[dof] - C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == 1
                    I[count] = ind; J[count] = ind; Val[count] = 2.0/dt + C2[dof] + C3[dof]; 
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count] + 4.0/dt # offset
                    center = count
                    count += 1
                    I[count] = ind; J[count] = indp; Val[count] = -C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == Nx[dof]
                    I[count] = ind; J[count] = indm; Val[count] = C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    I[count] = ind; J[count] = ind; Val[count] = 2.0/dt - C2[dof] + C3[dof]; 
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count] + 4.0/dt # offset
                    center = count
                    count += 1
                end
                # b
                p = b
                dof = 2
                indm = pos .- (0,1,0)
                indm = get_ind(Nx,indm)
                indp = pos .+ (0,1,0)
                indp = get_ind(Nx,indp)
                if 1 < p < Nx[dof] # inner
                    I[count] = ind; J[count] = indm; Val[count] = -C1[dof] + C2[dof] 
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    Val[center] += 2.0*C1[dof] + C3[dof] + 2.0/dt
                    Vc[center] = -Val[center] + 4.0/dt
                    I[count] = ind; J[count] = indp; Val[count] = -C1[dof] - C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == 1
                    Val[center] += 2.0/dt + C2[dof] + C3[dof];
                    Vc[center] = -Val[center] + 4.0/dt
                    I[count] = ind; J[count] = indp; Val[count] = -C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == Nx[dof]
                    I[count] = ind; J[count] = indm; Val[count] = C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    Val[center] = 2.0/dt - C2[dof] + C3[dof]; 
                    Vc[center] = -Val[center] + 4.0/dt
                end
                # c
                p = c
                dof = 3
                indm = pos .- (0,0,1)
                indm = get_ind(Nx,indm)
                indp = pos .+ (0,0,1)
                indp = get_ind(Nx,indp)
                if 1 < p < Nx[dof] # inner
                    I[count] = ind; J[count] = indm; Val[count] = -C1[dof] + C2[dof] 
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    Val[center] += 2.0*C1[dof] + C3[dof] + 2.0/dt
                    Vc[center] = -Val[center] + 4.0/dt
                    I[count] = ind; J[count] = indp; Val[count] = -C1[dof] - C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == 1
                    Val[center] += 2.0/dt + C2[dof] + C3[dof];
                    Vc[center] = -Val[center] + 4.0/dt
                    I[count] = ind; J[count] = indp; Val[count] = -C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                elseif p == Nx[dof]
                    I[count] = ind; J[count] = indm; Val[count] = C2[dof]
                    Ic[count] = I[count]; Jc[count] = J[count]; Vc[count] = -Val[count]
                    count += 1
                    Val[center] = 2.0/dt - C2[dof] + C3[dof]; 
                    Vc[center] = -Val[center] + 4.0/dt
                end
            end    
        end
    end    
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
    :l0 => -0.127061,
    :ve =>[0.2715618,0.0195524,0.0009720],
    :ka => [5.6772530,0.2520333,0.147],
    :si =>[0.0181427,0.0422960,0.034],
    :x0 => [0.05095958,0.06725220,0.00961570],
# harzard rate process (PSA params)
# ho(t) parameters
    :a => 0.024,
    :b => 1,
    # :b => 0.0000001, # test
    :T_asterisk => 2.5, # prepayment date,
    :gamma => 1,
    :k => 0.0007, # prepayment strike
    :k => 0.03, # prepayment strike
    # :k => 0.000000007, # prepayment strike
    :T_asterisk => 0.0, # prepayment date test
)

th = params[:ve]./params[:ka]
x_min = [0.0,0.0,0.0] # CIR
x_max = [0.1,0.1,0.1]
x_max = th*10.0
grid_d = Dict(
    :Nt => 30, # num of timesteps
    :Nx => [5,5,5], # num of state variables
    :Nx => [100,100,100], # num of state variables
    :x_min => x_min,
    :x_max => x_max,
)


l0 = params[:l0]
dx = (grid_d[:x_max] - grid_d[:x_min]) ./ (grid_d[:Nx] .- 1.0)
dt = params[:T] / grid_d[:Nt]

# Let's just write a script and think about how to organize it better


# Define grid. (Mind gamma is different depending on grid point)
Nx = grid_d[:Nx]
num_sv = length(Nx)
num_dof = prod(Nx)
num_nonzeros = get_n_nonzeros(Nx)
I = zeros(num_nonzeros)
J = zeros(num_nonzeros)
Val = zeros(num_nonzeros)
Ic = zeros(num_nonzeros)
Jc = zeros(num_nonzeros)
Vc = zeros(num_nonzeros)


get_mat_COO!(Nx, params, dx,dt,I,J,Val,Ic,Jc,Vc)
mat = sp.sparse(I,J,Val) # next step
matc = sp.sparse(Ic,Jc,Vc) # current step

T = params[:T]
tau_space = collect(0.0:dt:T)

# Try to solve it using an iterative IterativeSolvers

# 1. fQ(ru) = 1
fQ(ru) = 1.0
Q0 = zeros(num_dof) # initial value
Qc = zeros(num_dof) # current
Qn = zeros(num_dof) # next
for (ind,elem) in enumerate(Q0)
    Q0[ind] = fQ(elem)
end
compute_fd!(Q0,Qc,Qn,matc,mat,tau_space)
x0 = params[:x0]
x0_pos = round.((x0 -x_min)./dx) .+ 1
x0_ind = Int(get_ind(Nx,x0_pos))
println("Qn at t=0 ", Qn[x0_ind])
r0 = sum(x0) + l0


# 2. fH(ru) = h0(u) + gamma(k-ru)^+

# 3. fR(ru) = ru
fQ(ru) = ru
R0 = zeros(num_dof)
for (ind,elem) in enumerate(R0)
    pos = get_pos(Nx, ind)
    ru = sum(dx.*(pos.-1.0).+ x_min) + l0 
    R0[ind] = fQ(ru)
end