

include("mbs_functions.jl")


# include("/Users/jungjaeyong/projects/practices/quant/mbs/MBS_mc_cir.jl")
           
params = Dict(
# contract params
    :T => 30, # test
# Treasury
    :l0 => -0.127061,
    :l0 => 0.0, # test
    :ve => [0.2715618,0.0195524,0.0009720],
    :ka => [5.6772530,0.2520333,0.147],
    :si => [0.0181427,0.0422960,0.034],
    :x0 => [0.05095958,0.06725220,0.00961570],
    # :ve => [0.2715618,0.2715618,0.2715618], # test
    # :ka => [5.6772530,5.6772530,5.6772530], # test
    # :si => [0.0181427,0.0181427,0.0181427], # test
    # :x0 => [0.05095958,0.05095958,0.05095958], # test
# harzard rate process (PSA params)
# ho(t) parameters
    :a => 0.024,
    :b => 1,
    :b => 0.0, # test
    :T_asterisk => 2.5, # prepayment date,
    :gamma => 1,
    :gamma => 0, # test
    :k => 0.0007, # prepayment strike
    :k => 0.03, # prepayment strike
    :k => -1.0, # prepayment strike # test
    :T_asterisk => 1000.0, # prepayment date test # test
    )

ka = params[:ka]
si = params[:si]
ve = params[:ve]
th = ve./ka
s = 0
T = params[:T]
x0 = params[:x0]
num_states = length(x0)
num_paths = 3000
annual_steps = 1 # annual
# annual_steps = 12 # monthly
# annual_steps = 52 # weekly
num_steps = annual_steps * 30 
t_sim = collect(range(0,30,length = num_steps+1)) # initial at 0

# get intermediate variables for h(t)
ab = params[:a]*params[:b]
T_asterisk = params[:T_asterisk]
k = params[:k] # strike IR k
gamma = params[:gamma] # sensitivity to IR below strike k


println("start CIR sampling at T")
x_T = []
for i in 1:length(x0)
    push!(x_T, get_CIR_sample(num_paths,ka[i],si[i],th[i],s,T,x0[i]))
end

println("moments of x_T")
mean_x_T = stats.mean.(x_T)
theo_mean = get_CIR_mean.(T,ve,ka,x0)
mean_err = abs.(theo_mean - mean_x_T)./theo_mean *100.0
var_x_T = stats.var.(x_T)
std_x_T = sqrt.(var_x_T)
theo_var = get_CIR_var.(T,ve,ka,si,x0)
theo_std = sqrt.(theo_var)
std_err = abs.(std_x_T - theo_std)./theo_std *100.0
println("mean Error % ", mean_err)
println("std Error % ", std_err)

# get_CIR_sample_steps(num_paths,num_steps,ka[1],si[1],th[1],s,t,x0[1])

x_u = simulate_x(num_paths,num_steps,ka,si,th,s,T,x0) # x value at each step

# retrieve values
# you can add integration of shift also
r = get_r(x_u)
int_r = get_int_r(t_sim,x_u)
neg_int_r = -int_r # negative value is used

inter_h0_1 = zeros(num_steps+1) # same for every path
inter_h0_2 = zeros(size(r))
get_inter_h0!(inter_h0_1,inter_h0_2,t_sim,k,r,T_asterisk)
# inter_h0_1, inter_h0_2 = get_inter_h0(t_sim,k,r,T_asterisk)
int_h = zeros(num_paths,num_steps+1)

neg_int_r_h = zeros(size(neg_int_r))
@. neg_int_r_h = neg_int_r - int_h


# test whether those numbers makes sense
t,u = 0, T
R_0_T = get_R(t_sim,t,u,r,neg_int_r_h)
# compute Q
Q_0_T = get_Q(t_sim,t,u,neg_int_r_h)
# get_m


# Analytic Q
dl = ones(num_states)
ka = params[:ka]
si = params[:si]
tau = u-t
A = r_a_aux0(dl,th,ka,si,tau)
B =r_b_aux0(dl,ka,si,tau)
Q_analytic = exp(sum(A+B.*x0)) # without harzard
# I think get_Q is fine

# m equation
# numerator = int_0^T (1-exp(-m0*(T-u)))R(0,u)du
# denominator = int_0^T (1-exp(-m0*(T-u)))Q(0,u)du
m0 = 0.01 # initial guess


m = get_m_fi(m0,t_sim,T,r,neg_int_r_h)

# test parameters
# θ = 0.06, κ = 0.25, σ = 0.1
th = [0.06]
ka = [0.25]
si = [0.1]
x0 = [0.03] # 3,6,9%
k = x0[1] - 0.0025 # r - 0.0025 no shift

# rerun simulations
num_states = length(x0)
num_paths = 3000
# annual_steps = 1 # annual
# annual_steps = 12 # monthly
# annual_steps = 52 # weekly
num_steps = annual_steps * 30 
t_sim = collect(range(0,30,length = num_steps+1)) # initial at 0

x_u = simulate_x(num_paths,num_steps,ka,si,th,s,T,x0) # x value at each step

# you can add integration of shift also
r = get_r(x_u)
int_r = get_int_r(t_sim,x_u)
neg_int_r = -int_r # negative value is used

# no prepayment
a = 0.024
b = 0
ab = a*b
gamma = 0.0
T_asterisk = 1000.0
m0 = 0.01
get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
get_int_h!(int_h, ab, gamma, inter_h0_1, inter_h0_2)
@. neg_int_r_h = neg_int_r - int_h
# case 1
m_1 = zeros(size(t_sim))
for (idx,T) in enumerate(t_sim)
    if T <= 0.0
        continue
    end
    m_1[idx] = get_m_fi(m0,t_sim,T,r,neg_int_r_h)
end

# case 2 b=0.75, gamma = 0
a = 0.024
b = 0.75
ab = a*b
gamma = 0
T_asterisk = 2.5
get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
get_int_h!(int_h, ab, gamma, inter_h0_1, inter_h0_2)
@. neg_int_r_h = neg_int_r - int_h
m_2 = zeros(size(t_sim))
for (idx,T) in enumerate(t_sim)
    if T <= 0.0
        continue
    end
    m_2[idx] = get_m_fi(m0,t_sim,T,r,neg_int_r_h)
end


# case 3 b=0.0, gamma = 10
a = 0.024
b = 0.0
ab = a*b
gamma = 10.0
T_asterisk = 2.5
get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
get_int_h!(int_h, ab, gamma, inter_h0_1, inter_h0_2)
@. neg_int_r_h = neg_int_r - int_h
m_3 = zeros(size(t_sim))
for (idx,T) in enumerate(t_sim)
    if T <= 0.0
        continue
    end
    m_3[idx] = get_m_fi(m0,t_sim,T,r,neg_int_r_h)
end

# case 4 b=0.75, gamma = 10
a = 0.024
b = 0.75
ab = a*b
gamma = 10.0
T_asterisk = 2.5
get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
get_int_h!(int_h, ab, gamma, inter_h0_1, inter_h0_2)
@. neg_int_r_h = neg_int_r - int_h
m_4 = zeros(size(t_sim))
for (idx,T) in enumerate(t_sim)
    if T <= 0.0
        continue
    end
    m_4[idx] = get_m_fi(m0,t_sim,T,r,neg_int_r_h)
end
