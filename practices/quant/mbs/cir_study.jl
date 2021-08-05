
import Distributions; distr = Distributions
import Statistics; stats = Statistics
import NumericalIntegration; ni = NumericalIntegration
import Random

Random.seed!(3)


# include("/Users/jungjaeyong/projects/practices/quant/mbs/cir_study.jl")

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

function get_CIR_mean(t,ve,ka,x0)
    res = exp(-ka*t)*(x0 - ve/ka) + ve/ka
    return res
end

function get_CIR_var(t,ve,ka,si,x0)
    th = ve/ka
    res = si^2/ka*x0*(exp(-ka*t)-exp(-2.0*ka*t)) +th*si^2/(2.0*ka)*(1-exp(-ka*t))^2
    return res
end

function get_CIR_sample!(buffer,num_paths,ka,si,th,s,t,x0)
    kappa = ka
    gamma = si
    vbar = th
    v_s = x0 # initial state variable
    delta = 4.0 *kappa*vbar/gamma^2
    c= 1.0/(4.0*kappa)*gamma^2*(1.0-exp(-kappa*(t-s)))
    kappaBar = 4.0*kappa*v_s*exp(-kappa*(t-s))/(gamma^2*(1.0-exp(-kappa*(t-s))))
    d = distr.NoncentralChisq(delta,kappaBar)
    # @show d,buffer,rand(d, num_paths)
    buffer .= rand(d, num_paths)
    buffer .*= c
    return buffer
end

function get_CIR_sample(num_paths,ka,si,th,s,t,x0::Float64)
    if num_paths == 1
        buffer = [NaN]
        get_CIR_sample!(buffer,num_paths,ka,si,th,s,t,x0)
        return buffer[1]
    else
        buffer = zeros(num_paths)
        get_CIR_sample!(buffer,num_paths,ka,si,th,s,t,x0)
        return buffer
    end
end


function get_CIR_sample(num_paths,ka,si,th,s,t,x0::Vector)
    buffer = zeros(num_paths)
    get_CIR_sample!(buffer,num_paths,ka,si,th,s,t,x0)
    return buffer
end

function get_CIR_sample_steps(num_paths,num_steps,ka,si,th,s,t,x0)
    dt = (t-s) / num_steps
    X = zeros(num_paths, num_steps+1)
    X[:,1] .= x0 # step-1 is initial

    for i in 1:num_steps
        for j in 1:num_paths
            X[j,i+1] = get_CIR_sample(1,ka,si,th,0.0,dt,X[j,i])
        end
    end
    return X
end

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
t = params[:T]
x0 = params[:x0]
num_paths = 2000
annual_steps = 12 # monthly
annual_steps = 52 # weekly
num_steps = annual_steps * 30 
t_sim = collect(range(0,30,length = num_steps+1))

println("start CIR sampling at T")
x_T = []
for i in 1:length(x0)
    push!(x_T, get_CIR_sample(num_paths,ka[i],si[i],th[i],s,t,x0[i]))
end

println("moments of x_T")
mean_x_T = stats.mean.(x_T)
theo_mean = get_CIR_mean.(t,ve,ka,x0)
mean_err = abs.(theo_mean - mean_x_T)./theo_mean *100.0
var_x_T = stats.var.(x_T)
theo_var = get_CIR_var.(t,ve,ka,si,x0)
var_err = abs.(theo_var - var_x_T)./theo_var *100.0
println("mean Error % ", mean_err)
println("var Error % ", var_err)

# get_CIR_sample_steps(num_paths,num_steps,ka[1],si[1],th[1],s,t,x0[1])

x_u = [] # x value at each step
for i in 1:length(x0)
    push!(x_u, get_CIR_sample_steps(num_paths,num_steps,ka[i],si[i],th[i],s,t,x0[i]))
    # push!(x_u, get_CIR_sample_steps(num_paths,num_steps,ka[i],si[i],th[i],s,t,x0[i])[:,2:end])
end

@time begin
ni.integrate(t_sim, @view x_u[1][1,:])
end
    