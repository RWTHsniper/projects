import Distributions; dist = Distributions
import Random
import Statistics; stats = Statistics
import NumericalIntegration; ni = NumericalIntegration
import Optim

# include("C://Users//golde//Documents//GitHub//projects//practices//quant//mbs//mc.jl")

# Price of MBS
# if t< tau < T
# Mt = discounted c + discounted P_tau 
# if T <= tau
# Mt = 

function get_CIR_mean(t,ve,ka,si,x0)
    res = exp(-ka*t)*(x0 - ve/ka) + ve/ka
    return res
end

function get_P_no_prepay(P0,m,t,T)
# Principle without prepayment
    Pt = P0*(1-exp(-m*(T-t)))/(1-exp(-m*T))
    return Pt
end

function get_c(m,P0,T)
    # without prepayment
    res = m*P0/(1.0-exp(-m*T))
    return res
end

function get_ho_t(a,b,t,T_asterisk)
# ho(t) deterministic
    res = b*a*max(t,T_asterisk)
    return res
end

function get_h_t(a,b,gamma,k,r,t,T_asterisk)
    ho_t = get_ho_t(a,b,t,T_asterisk)
    res = ho_t + gamma*max(k-r,0.0)
    return res
end

function get_h_t(a,b,gamma,k,r::Union{Vector,SubArray},t,T_asterisk)
    ho_t = get_ho_t(a,b,t,T_asterisk)
    res = ho_t .+ gamma*max.(k.-r,0.0) # max(k-r,0)
    return res
end

get_P_no_prepay(1.0,0.03,11,10)

# Model parameters
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
e_dist = dist.Exponential(1) # exponential distribution about prepayment. Integration of harzard process > e

# Let's write a simulation code
num_paths= 5
annual_steps = 12 # steps per year
num_paths= 1000
annual_steps = 300 # steps per year
dt = 1/annual_steps
T = 10 # year
# T = 30 # test
num_steps = T * annual_steps

# variables in simulations
num_x = length(params[:x0])
x = zeros(num_paths,num_x,num_steps)
r = zeros(num_paths,num_steps) # r is sum of x
h = zeros(num_paths,num_steps)
e = zeros(num_paths,num_steps) # random threshold for prepayment
Random.rand!(e_dist,e) # generate e before running simulation since it is simpler and independent
d_B = dist.Normal(0.0,1.0) # standard distribution
dW_x = zeros(num_paths,num_x) # Brownian motion for x

# Initialization
Random.seed!(3)

ve = params[:ve]
ka = params[:ka]
si = params[:si]
l0 = params[:l0]
a = params[:a]
b = params[:b]
k = params[:k]
gamma = params[:gamma]
T_asterisk = params[:T_asterisk]

current = 1 # current step
# compute initial interest rate
r[:,current] .= params[:l0]
for i in 1:num_x
    x[:,i,current] .= params[:x0][i] # assign initial treasury variable
    r[:,current] .+= params[:x0][i] # initial interest rate
end
println("Initial r ", r[1,1])
h[:,current] = get_h_t(a,b,gamma,k,r[:,current],0.0,T_asterisk)
t_sim = collect(0:num_steps-1)*dt

# run simulation

# Run simulation at each timestep. compute for the next value
for current in 1:(num_steps-1)
    global dW_x
    t = dt*current
    # simulate treasury parameters
    Random.rand!(d_B,dW_x) # Brownian motions for x
    mean_dW_x = stats.mean(dW_x,dims=1)
    std_dW_x = stats.std(dW_x,dims=1)
    dW_x .-= mean_dW_x
    dW_x ./= std_dW_x
    dW_x *= sqrt(dt) # std of Brownian motion
    next = current + 1
    for i in 1:num_x
        drift = (ve[i] .- ka[i]*x[:,i,current])*dt
        vol = si[i]*sqrt.(x[:,i,current]).*dW_x[:,i]
        vol += 0.25*si[i]^2*(dW_x[:,i].^2 .- dt) # Milstein scheme
        x[:,i,next] = x[:,i,current] + drift + vol
        x[:,i,next] = max.(0.0, x[:,i,next]) # reflection for CIR process
    end
    
    # compute r
    r[:,next] .= l0
    for i in 1:num_x
        r[:,next] .+= x[:,i,next] # initial interest rate
    end
    
    # compute harzard rate
    global h_t
    x_view = @view x[:,end,next]
    h[:,next] = get_h_t(a,b,gamma,k,x_view,t,T_asterisk)
end

println("Check whether CIR process simulation is reliable")
# statistics for x
mean_x = stats.mean(x[:,:,end],dims=1)[:]
theo_mean_x = get_CIR_mean.(T,ve,ka,si,params[:x0])
@show mean_x
@show theo_mean_x

# Statistics for r
mean_r = stats.mean(r[:,end])
std_r = stats.std(r[:,end])
theo_mean_r = sum(theo_mean_x) + params[:l0]
@show mean_r, std_r
@show theo_mean_r

# integrate h_t
int_h = zeros(size(h))
# Trapezoid scheme
for j in 1:num_steps-1
    int_h[:,j+1] = int_h[:,j] + dt/2.0*(h[:,j]+h[:,j+1])
end
prepay_flag = int_h .> e
prepay_step = zeros(Int64, num_paths)
for i in 1:length(prepay_step)
    try
        prepay_step[i] = findfirst(prepay_flag[i,:])
    catch e
        prepay_step[i] = num_steps # prepayment did not happen
    end
end

println("maximum prepayment step ", maximum(prepay_step), " total steps ", num_steps)
println("minimum prepayment step ", minimum(prepay_step), " total steps ", num_steps)

# First, we need to compute m(0,T)
P0 = 1.0 # presume initial mortgage price
# m = 0.03 # mortgage rate


# function get_int_t_u_r(t_sim,r,t_ind,u_ind) 
#     t = t_sim[t_ind]
#     u = t_sim[u_ind]
#     res = ni.integrate(t_sim[t_ind:u_ind], r[t_ind:u_ind])
# end

int_0_u_r = zeros(num_paths,num_steps) # int_0^(u) rs ds at each step
int_0_u_h = zeros(num_paths,num_steps) # int_0^(u) hs ds at each step
for i in 1:num_paths
    for j in 2:num_steps
        int_0_u_r[i,j] = ni.integrate(view(t_sim,1:j), view(r,i,1:j))
        int_0_u_h[i,j] = ni.integrate(view(t_sim,1:j), view(h,i,1:j))
    end
end
# for example, if you want to compute int_u1_u2_r, then compute int_0_u2_r - int_0_u1_r
# There is no need to compute numerical integration over r and h anymore

int_0_u_r_h = int_0_u_r + int_0_u_h
function get_Q(int_t_u_r_h) # (3.2)
    # get int_r_h_u at u
    res = exp.(-int_t_u_r_h)
    res = stats.mean(res)
    return res
end
# get_Q(int_t_u_r_h)

function get_R(int_t_u_r_h, r_u) # 3.7
    res = r_u.*exp.(-int_t_u_r_h)
    res = stats.mean(res)
    return res
end
# get_R(int_t_u_r_h, r_u)

# calculate m
m0 = 0.01 # initial guess on m
function m_func(m0,t_sim,T,int_0_u_r_h,r) # 3.9
    d_num = zeros(size(t_sim))
    d_den = zeros(size(t_sim))
    for (j,u) in enumerate(t_sim)
        r_u = view(r,:,j)
        factor = (1-exp(-m0*(T-u)))
        d_num[j] = factor*get_R(int_0_u_r_h[:,j], r_u)
        d_den[j] = factor*get_Q(int_0_u_r_h[:,j])
    end
    num = ni.integrate(t_sim, d_num)
    den = ni.integrate(t_sim, d_den)
    m = num/den
    return m
end

obj_fun = x -> (x - m_func(x,t_sim,T,int_0_u_r_h,r))^2
res = Optim.optimize(obj_fun, 0.001,0.04)
m = res.minimizer
println("Mortgage rate ",m)
c = get_c(m,P0,T) # coupon (2.1)

# Computation of M0 (3.1)
M0 = zeros(num_paths)
buffer = zeros(num_steps) # buffer to save intermediate result
for i in 1:num_paths # looping over each path
    p_step = prepay_step[i] # prepayment step for i-th path
    t_view = view(t_sim, 1:p_step) # use view to avoid deepcopy
    int_r = int_0_u_r[i,p_step]
    if p_step < num_steps # prepayment occured
        discount = exp(-int_r)
        M0[i] = discount*get_P_no_prepay(P0,m,t_sim[p_step],T)
    end
    # compute cashflow from coupon
    # compute discounted c
    for (j, t) in enumerate(t_view) # j-time step
        discount = exp(-int_r)
        buffer[j] = discount
    end
    M0[i] += c * ni.integrate(t_view, view(buffer,1:p_step)) # mind deep copy
end

println("M0 at t=0 ",stats.mean(M0))

# Another formulation for M




# Test for M0
# The following expression should be zero 3.6
t = 0.0
t_ind = 1
# buffer = zeros(length(t_sim)-t_ind+1)
buffer = zeros(length(t_sim))
for (u_ind,u) in enumerate(t_sim)
    factor = (1.0-exp(-m*(T-u)))/(1.0-exp(-m*(T-t)))
    int_view = view(int_0_u_r_h,:,t_ind:u_ind)
    r_view = view(r,:,u_ind)
    factor2 = m * get_Q(int_view) - get_R(int_view, r_view)
    buffer[u_ind] = factor*factor2
end
int_test = ni.integrate(t_sim, buffer)
println("int_test is supposed to be 0 at t=0", int_test)