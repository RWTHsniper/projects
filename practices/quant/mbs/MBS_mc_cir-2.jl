

include("mbs_functions.jl")

import Plots; plt = Plots
import Optim
import PyPlot

# include("/Users/jungjaeyong/projects/practices/quant/mbs/MBS_mc_cir-2.jl")
           
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
    # MC params
    :annual_steps => 1, # annual
    :num_paths => 3000,
    )

# test parameters
params[:th] = [0.06]
params[:ka] = [0.25]
params[:si] = [0.1]
params[:x0] = [0.03] # 3,6,9%
params[:k] = params[:x0][1] - 0.0025 # r - 0.0025 no shift
params[:a] = 0.024

mutable struct MBS_cir
    # CIR
    l0::Float64
    dl::Vector{Float64}
    ve::Vector{Float64}
    ka::Vector{Float64}
    si::Vector{Float64}
    x0::Vector{Float64}
    th::Vector{Float64}
    num_states::Real
    # contract
    t0::Real # initial time
    T::Real # maturity
    # harzard process
    a::Real
    b::Real
    T_asterisk::Real
    gamma::Real
    k::Real
    # MC parameters and variables
    annual_steps::Real
    num_paths::Real
    num_steps::Real
    t_sim::Vector{Float64}
    x::Vector # CIR sampling at each steps
    r
    neg_int_r
    inter_h0_1 # intermediate value
    inter_h_r # intermediate value
    int_h
    neg_int_r_h
    function MBS_cir(params::Dict)
        this = new()
        for (key,val) in params
            # @show key,val
            setfield!(this,key,val)
        end
        this.num_states = length(this.ve)
        this.dl = ones(this.num_states)
        if !(:t0 in keys(params))
            this.t0 = 0.0
        end
        if !(:th in keys(params))
            this.th = this.ve ./ this.ka
        end
        num_steps = this.annual_steps*this.T
        num_paths = this.num_paths
        this.num_steps = num_steps
        this.t_sim = collect(range(0,this.T,length = num_steps+1)) # initial at 0
        this.r = zeros(num_paths,num_steps+1)
        this.neg_int_r = zeros(num_paths,num_steps+1)
        this.inter_h0_1 = zeros(num_steps+1)
        this.inter_h_r = zeros(num_paths,num_steps+1)
        this.int_h = zeros(num_paths,num_steps+1)
        this.neg_int_r_h = zeros(num_paths,num_steps+1)
        return this
    end
end

function simulate_x!(m::MBS_cir, shift=nothing)
    m.x = simulate_x(m.num_paths,m.num_steps,m.ka,m.si,m.th,m.t0,m.T,m.x0)
    m.r = get_r(m.x)
    m.neg_int_r = - get_int_r(m.t_sim,m.x)
    return 
end

function compute!(m::MBS_cir,calibration=false)
    # compute intermediate variables
    if calibration == false # unnecessary for prepayment intensity calibration
        get_inter_h0_1!(m.inter_h0_1,m.t_sim,m.T_asterisk)
    end
    get_inter_h_r!(m.inter_h_r,m.t_sim,m.k,m.r)
    ab = m.a*m.b
    get_int_h!(m.int_h, ab, m.gamma, m.inter_h0_1, m.inter_h_r)
    @. m.neg_int_r_h = m.neg_int_r - m.int_h
end

function get_m_curve(m::MBS_cir)
    m_curve = zeros(size(m.t_sim))
    m0 = sum(m.th) # steady-state value as an initial guess
    # m0 = 0.01
    for (idx,s) in enumerate(m.t_sim)
        if s <= 0.0
            continue
        end
        m_curve[idx] = get_m_fi(m0,m.t_sim,s,m.r,m.neg_int_r_h) # fixed-iteration
    end
    return m_curve
end

function update!(m::MBS_cir,calibration=false;kwargs...)
    # update! does not re-simulate the paths
    for (key,val) in kwargs
        setfield!(m,key,val)
    end
    compute!(m::MBS_cir,calibration)
end

mbs_model = MBS_cir(params)
simulate_x!(mbs_model)
compute!(mbs_model)

cases = [Dict(:a=>0.024, :b=>0.0, :gamma=>0.0, :T_asterisk=>1000.0), # no prepayment
        Dict(:a=>0.024, :b=>0.75, :gamma=>0.0, :T_asterisk=>2.5),
        Dict(:a=>0.024, :b=>0.0, :gamma=>10.0, :T_asterisk=>2.5),
        Dict(:a=>0.024, :b=>0.75, :gamma=>10.0, :T_asterisk=>2.5)] 

m_curves = []
for i in 1:length(cases)
    update!(mbs_model; cases[i]...)
    push!(m_curves,get_m_curve(mbs_model))
end

# Use PyPlot instead of Plots
PyPlot.plot(mbs_model.t_sim[2:end],m_curves[1][2:end],label="No prep.")
for i in 2:length(cases)
    PyPlot.plot(mbs_model.t_sim[2:end],m_curves[i][2:end],label="b="*string(cases[i][:b])*" gamma="*string(cases[i][:gamma]))
end
PyPlot.title("r_0 = "*string(mbs_model.x0[1]))
PyPlot.legend()

# p = plt.plot(mbs_model.t_sim, m_curves[1])
# plt.plot!(p, mbs_model.t_sim, m_curves[2])
# plt.plot!(p, mbs_model.t_sim, m_curves[3])
# plt.plot!(p, mbs_model.t_sim, m_curves[4])
# for i in 1:length(cases)
#     if i==1
#         p = plt.plot(mbs_model.t_sim, m_curves[i])
#         global p
#     else
#         plt.plot!(p, mbs_model.t_sim, m_curves[i])
#     end
# end

# Fit parameters to a mortgage-rate curve

"""
# Arguments
- `x`: vector of variables (b and gamma)
- `m`: MBS model
- `mat`: maturities from the market
- `m_market`: mortgage rate curve from the market
"""
function objective_function(x, m::MBS_cir, mat,m_market)
    update!(m, true; b=x[1], gamma=x[2])
    m0 = 0.01;
    res = 0.0
    for (idx,T) in enumerate(mat)
        # @show m.t_sim, T
        if T <= 0.0
            continue
        end
        m_idx = get_m_fi(m0,m.t_sim,T,m.r,m.neg_int_r_h)
        res += (m_idx - m_market[idx])^2
    end
    if x[1] < 0.0 # b
        res += abs(x[1])*100 # penalty
    end
    if x[2] < 0.0 # gamma 
        res += abs(x[1])*100 # penalty
    end
    return res
end

initial_guess = [mbs_model.b, mbs_model.gamma]
idx = 3
obj = x -> objective_function(x,mbs_model,mbs_model.t_sim,m_curves[idx]) # let's try to fit curve idx
res = Optim.optimize(obj, initial_guess)
sol = res.minimizer
mbs_model.b = sol[1]; mbs_model.gamma = sol[2]
calibrated_curve = get_m_curve(mbs_model)

PyPlot.plot(mbs_model.t_sim[2:end],m_curves[idx][2:end],label=string(idx)*" curve")
PyPlot.plot(mbs_model.t_sim[2:end],calibrated_curve[2:end],label="calibrated curve")
PyPlot.legend()
