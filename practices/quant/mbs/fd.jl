
# include("/Users/jungjaeyong/projects/practices/quant/mbs/fd.jl")


import SparseArrays; sp = SparseArrays
import IterativeSolvers; is = IterativeSolvers



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

x_min = [0.0,0.0,0.0] # CIR
x_max = [0.1,0.1,0.1]
grid_d = Dict(
    :Nt => 30, # num of timesteps
    :Nx => [5,5,5], # num of state variables
    :x_min => x_min,
    :x_max => x_max,
)



dx = (grid_d[:x_max] - grid_d[:x_min]) ./ (grid_d[:Nx] .- 1.0)


# Let's just write a script and think about how to organize it better


# Define grid. (Mind gamma is different depending on grid point)
Nx = grid_d[:Nx]
num_sv = length(Nx)
num_nonzeros = num_sv*2 + sum(Nx)
I = zeros(num_nonzeros)
J = zeros(num_nonzeros)
Val = zeros(num_nonzeros)

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

counter = 1
for c in 1:Nx[3]
    for b in 1:Nx[2]
        for a in 1:Nx[1]
        global counter,Nx
            # For a grid (a,b,c), we can build up I,J,Val
            pos = (a,b,c)
            ind = get_ind(Nx, pos)
            # Let's think about how to efficiently construct an algorithm
            if a == 1
                ip1 = get_ind(2,b,c) # ip1
            else if a == Nx[1]
            end
            # construct row by row
            @show ind
            counter += 1
        end
    end
end
@show counter

# Let's go


# Define initial condition (tau=0, t=T)


# Try to solve it using an iterative IterativeSolvers

