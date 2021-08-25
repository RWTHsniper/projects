
#=
include("/Users/jungjaeyong/projects/practices/quant/mbs/fd_functions.jl")
=#

import SparseArrays; sp = SparseArrays
import IterativeSolvers; is = IterativeSolvers
import LinearAlgebra; la = LinearAlgebra
import Interpolations; ip = Interpolations

include("cir_functions.jl")

"""
Takes real and vector inputs
"""
function get_x(x_min,dx,ind)
    return x_min + dx.*(ind.-1)
end

"""
Takes real and vector inputs
"""
function get_x_ind(x,x_min,dx)
    ind = round.(Int,(x-x_min)./dx .+ 1)
    return ind
end

function get_fxx(dx::Real,order=2)
    if order == 2 # f(j-1),f(j),f(j+1)
        res = [1.0, -2.0, 1.0] / dx^2
    else
        res = nothing
    end
    return res
end

function get_fx(dx::Real,order=2)
    if order == 2
        res = [-1.0, 0.0, 1.0] / (2.0*dx) # f(j-1),f(j),f(j+1)
    elseif order == 1
        res = [-1.0, 1.0] / (dx) # f(j-1),f(j) or f(j), f(j+1)
    end
    return res
end

"""
Computes the RHS for a second-order PDE with zero-gamma boundary conditions.
"""
function get_mat_discretized!(mat::la.Tridiagonal,x_min::Real,dx::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    mat .= 0.0
    num_rows = size(mat,1)
    for ind in 1:num_rows # row index
        x = get_x(x_min,dx,ind)
        indm = ind - 1
        indp = ind + 1
        if !isnothing(get_coeff_0)
            mat[ind,ind] += get_coeff_0(t,x) # always we have it discount term
        end
        if ind == 1 # bc
            # zero gamma: coeff_2 = 0
            # fx
            if !isnothing(get_coeff_1)
                mat[ind,ind:indp] .+= get_coeff_1(t,x) * get_fx(dx,1)
            end
        elseif ind == num_rows # bc
            # zero gamma: coeff_2 = 0
            # fx
            if !isnothing(get_coeff_1)
                mat[ind,indm:ind] .+= get_coeff_1(t,x) * get_fx(dx,1)
            end
        else # inner node
            if !isnothing(get_coeff_2)
                mat[ind,indm:indp] .+= get_coeff_2(t,x) * get_fxx(dx)
            end
            if !isnothing(get_coeff_1)
                mat[ind,indm:indp] .+= get_coeff_1(t,x) * get_fx(dx) 
            end
        end
    end
end


function get_mat_explicit!(mat_explicit::la.Tridiagonal,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    mat_explicit .= 0.0
    get_mat_discretized!(mat_explicit,x_min,dx,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    for ind in 1:size(mat_explicit,1)
        mat_explicit[ind,ind] += 1.0/dt
    end
end

function get_mat_explicit(n_dof,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    mat_explicit = la.Tridiagonal(zeros(n_dof-1),zeros(n_dof),zeros(n_dof-1)) # matrix for nodes at current step
    get_mat_explicit!(mat_explicit,x_min,dx,dt,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    return mat_explicit
end

function get_mat_implicit!(mat_implicit::la.Tridiagonal,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    get_mat_discretized!(mat_implicit,x_min,dx,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    mat_implicit .= -mat_implicit
    for ind in 1:size(mat_implicit,1) # row index
        mat_implicit[ind,ind] += 1.0/dt
    end
end

function get_mat_implicit(n_dof,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    mat_implicit = la.Tridiagonal(zeros(n_dof-1),zeros(n_dof),zeros(n_dof-1)) # matrix for nodes at current step
    get_mat_implicit!(mat_implicit,x_min,dx,dt,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2) # reuse code
    return mat_implicit
end

function gat_mat_CN!(mat_next,mat_current,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    get_mat_discretized!(mat_current,x_min,dx,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    mat_current .*= 0.5 
    mat_next .= -mat_current
    t_diag = 1.0/dt # diagonal term for time discretization
    for ind in 1:size(mat_current,1) # row index
        mat_current[ind,ind] += t_diag
        mat_next[ind,ind] += t_diag
    end
end

function gat_mat_CN(n_dof,x_min::Real,dx::Real,dt::Real,t=nothing;get_coeff_0=nothing,get_coeff_1=nothing,get_coeff_2=nothing)
    mat_next = la.Tridiagonal(zeros(n_dof-1),zeros(n_dof),zeros(n_dof-1)) # matrix LHS
    mat_current = la.Tridiagonal(zeros(n_dof-1),zeros(n_dof),zeros(n_dof-1)) # matrix for RHS
    gat_mat_CN!(mat_next, mat_current,x_min,dx,dt,t;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    return mat_next, mat_current
end

#=
Supported discretization schemes

1. Explicit
2. Implicit
3. Crank-Nicolson
4. ADI (Alternating-direction implicit method)

Todo
1. Scheme initializer (mem allocation)
2. Discretization calculator

=#
function main()
    params = Dict(
    # contract params
        :T => 30, # test
    # Treasury
        :l0 => 0.0,
        :ve => 0.0027646771594520936,
        :ka => 0.06491552810007564, 
        :si => 0.03362592979733442, 
        :x0 => 0.022, # test
    # harzard rate process (PSA params)
    # ho(t) parameters
        :a => 0.0,
        :b => 0.0, 
        :gamma => 20.0,
        :gamma => 0.0,
        :k => 0.02, # prepayment strike
        # :k => 0.015, # prepayment strike
        :T_asterisk => 1000.0, # prepayment date test # test
        )

    grid_d = Dict(
        :Nt => params[:T]*12, # num of timesteps
        :Nt => params[:T]*12*5*7, # num of timesteps
        # :Nx => [3,3,3], # num of state variables
        :Nx => 8, # num of state variables
        # :Nx => 64, # num of state variables
        :Nx => 256, # num of state variables 
        # :Nx => 512, # num of state variables. Dont work for explicit
        # :Nx => 1024, # num of state variables. Dont work for explicit
        # :Nx => [128,128,32], # num of state variables
        :x_min => 0.0000001,
        :x_max => 0.6,
    )
        

    ka = params[:ka]
    si = params[:si]
    ve = params[:ve]
    x0 = params[:x0]
    gamma = params[:gamma]
    k = params[:k]

    dx = (grid_d[:x_max] - grid_d[:x_min]) / (grid_d[:Nx] - 1.0)
    dt = params[:T] / grid_d[:Nt]
    Nx = grid_d[:Nx]
    x_min = grid_d[:x_min]


    # compose the explicit scheme
    n_dof = Nx[1]

    # coefficients depending on a model
    get_coeff_2(t,x) = 0.5*si^2*x # 2nd order coefficient
    get_coeff_1(t,x) = ve - ka*x
    function get_coeff_0(t,x)
        if x <= k
            res = -(x + gamma*(k-x))
        else
            res = -x
        end
        return res
    end

    # matrices for Crank-Nicolson scheme


    T = params[:T]
    tau_space = collect(0.0:dt:T)

    # Try to solve it using an iterative IterativeSolvers

    # 1. fQ(ru) = 1
    fQ(ru) = 1.0
    Q0 = zeros(Nx[1]) # initial value
    Qn = zeros(Nx[1]) # next
    for (ind,elem) in enumerate(Q0)
        Q0[ind] = fQ(elem)
    end

    mat_explicit = get_mat_explicit(n_dof,x_min,dx,dt;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    Qc = deepcopy(Q0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Qn .= dt*mat_explicit*Qc
        Qc .= Qn
    end

    # Q value
    x0_ind = get_x_ind(x0,x_min,dx)
    println("Explicit scheme: Qn at t=0 ", Qn[x0_ind]) # It is quite sensitive to discretizaitons

    # Implicit scheme
    mat_implicit = get_mat_implicit(n_dof,x_min,dx,dt;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    Qc = deepcopy(Q0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Qn .= mat_implicit\(Qc/dt)
        Qc .= Qn
    end
    # Q value
    x0_ind = get_x_ind(x0,x_min,dx)
    println("Implicit scheme: Qn at t=0 ", Qn[x0_ind])

    # Crank-Nicolson scheme.
    mat_next, mat_current = gat_mat_CN(n_dof,x_min,dx,dt;get_coeff_0=get_coeff_0,get_coeff_1=get_coeff_1,get_coeff_2=get_coeff_2)
    Qc = deepcopy(Q0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Qn .= (mat_next)\(mat_current*Qc) # explicit
        Qc .= Qn
    end
    # Q value
    x0_ind = get_x_ind(x0,x_min,dx)
    println("CN scheme: Qn at t=0 ", Qn[x0_ind])


    # Analytic Q
    dl = 1
    ka = params[:ka]
    si = params[:si]
    th = params[:ve] ./ params[:ka]
    tau = T
    A = r_a_aux0(dl,th,ka,si,tau)
    B =r_b_aux0(dl,ka,si,tau)
    bond_price = exp(sum(A+B*x0))
    println("Analytic solution for bond price: ",bond_price)


    # R
    fR(ru) = ru
    R0 = zeros(Nx) # initial value
    Rn = zeros(Nx) # next
    for (ind,elem) in enumerate(R0)
        x_ind = get_x(x_min,dx,ind)
        R0[ind] = fR(x_ind)
    end

    # explicit
    Rc = deepcopy(R0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Rn .= dt*mat_explicit*Rc
        Rc .= Rn
    end
    x0_ind = get_x_ind(x0,x_min,dx)
    println("Explicit scheme: Rn at t=0 ", Rn[x0_ind])

    # implicit
    Rc = deepcopy(R0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Rn .= mat_implicit\(Rc/dt)
        Rc .= Rn
    end
    # Q value
    x0_ind = get_x_ind(x0,x_min,dx)
    println("Implicit scheme: Rn at t=0 ", Rn[x0_ind])

    # CN scheme
    Rc = deepcopy(R0) # current
    for (ind,tau) in enumerate(tau_space)
        if ind == length(tau_space)
            continue
        end
        Rn .= (mat_next)\(mat_current*Rc) # explicit
        Rc .= Rn
    end
    # Q value
    x0_ind = get_x_ind(x0,x_min,dx)
    println("CN scheme: Rn at t=0 ", Rn[x0_ind])

    # Analytic R
    dl = 1.0
    ka = params[:ka]
    si = params[:si]
    th = params[:ve] ./ params[:ka]
    tau = T
    A = r_a_aux0(dl,th,ka,si,tau)
    A_du = r_a_aux0_du(dl,th,ka,si,tau)
    B =r_b_aux0(dl,ka,si,tau)
    B_du = r_b_aux0_du(dl,ka,si,tau)
    bond_price = exp(sum(A+B*x0))
    R_anal = (A_du + B_du*x0)*bond_price
    println("Analytic solution for R: ", R_anal)

end

# main()