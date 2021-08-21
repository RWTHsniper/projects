

include("mbs_functions.jl")

import Plots; plt = Plots
import Optim
#import PyPlot # Not using at Conning
import Interpolations; ip = Interpolations

# include("/Users/jungjaeyong/projects/practices/quant/mbs/MBS_mc_cir-2.jl")
#include("C:\\Users\\Jaejun\\Documents\\conning\\MBS\\codes\\jj\\MBS_mc_cir-2.jl")

# kill a,b, shift which are deterministic parts
tot_params = Dict(
# contract params
    :T => 30, # test
# Treasury
    :l0 => -0.127061,
    :l0 => 0.0, # test
    :ve => [0.0027646771594520936], # test
    :ka => [0.06491552810007564], # test
    :si => [0.03362592979733442], # test
    :x0 => [0.022], # test
# harzard rate process (PSA params)
# ho(t) parameters
    :a => 0.0,
    :b => 0.0, # test 
    :gamma => 20.0,
    :k => 0.02, # prepayment strike
    :k => 0.001, # prepayment strike test
    :T_asterisk => 1000.0, # prepayment date test # test
    # MC params
    :annual_steps => 1, # annual
    :annual_steps => 12, # monthly
    # :annual_steps => 52, # weekly
    :num_paths => 10000,
    )

# we are going to only work on the last term
params = copy(tot_params)
last_ind = [length(tot_params[:x0])]
params[:ve] = tot_params[:ve][last_ind]
params[:ka] = tot_params[:ka][last_ind]
params[:si] = tot_params[:si][last_ind]
params[:x0] = tot_params[:x0][last_ind]

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
    shift_f
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
    int_shift
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
        this.int_shift = zeros(num_steps+1)
        get_int_shift!(this.int_shift,this.shift_f,this.t_sim)
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
    # get_inter_h_r!(m.inter_h_r,m.t_sim,m.k,m.r)
    get_inter_h_r!(m.inter_h_r,m.t_sim,m.k,m.x[end]) # use m.x[end] for the harzard process
    ab = m.a*m.b
    get_int_h!(m.int_h, ab, m.gamma, m.inter_h0_1, m.inter_h_r)
    for j in 1:m.num_steps
        for i in 1:m.num_paths
            m.neg_int_r_h[i,j] = m.neg_int_r[i,j] - m.int_h[i,j] - m.int_shift[j]
        end
    end
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
        if key == :shift_f
            get_int_shift!(m.int_shift,m.shift_f,m.t_sim) # update integration of shift
        end
    end
    compute!(m::MBS_cir,calibration)
end

# shift = [0.0 -0.01; 0.25 -0.01; 0.5 -0.01; 1.0 -0.01; 30.0 -0.01]
shift = [0 -0.2710956205188755; 0.25 -0.2741336132511397; 0.5 -0.27323142762050856; 1 -0.27830959346297984; 2 -0.28042256589576214; 3 -0.2811959066711729; 4 -0.28777025237033566; 5 -0.28422889483148406; 6 -0.28601809510984444; 7 -0.2913259018579911; 8 -0.28263372606069165; 9 -0.2953097574503388; 10 -0.2890850622651007; 11 -0.28777512679538836; 12 -0.293330906716708; 13 -0.296742825348417; 14 -0.2903551259059402; 15 -0.29107759881067663; 16 -0.30104630234126345; 17 -0.2923522127025746; 18 -0.3043166208707728; 19 -0.2846992095304273; 20 -0.3128582063936078; 21 -0.292034817101383; 22 -0.2987534978293915; 23 -0.3036266481509428; 24 -0.29544552575531624; 25 -0.3031655691557211; 26 -0.2996322925586268; 27 -0.302332700261461; 28 -0.30225313316559954; 29 -0.30079838595288055; 30 -0.30256867900705403]
shift[:,2] .= 0.0 # test
params[:shift_f] = ip.LinearInterpolation(shift[:,1],shift[:,2])

mbs_model = MBS_cir(params)
simulate_x!(mbs_model)
compute!(mbs_model)


annualtimes = vcat([0.25,0.5,1,2,3,4,5,7,8,9,10,15,20,25,30]) # times at which you want to calculate the nortgage rate Q ,R for the ouptput
tm = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
Q_mc = get_Q(mbs_model.t_sim,mbs_model.t_sim,mbs_model.T,mbs_model.neg_int_r_h)
Q_f = ip.LinearInterpolation(mbs_model.t_sim,Q_mc) # Q function
Q_tm = Q_f(tm) # Q value at tm

m = mbs_model
zcb_price = get_zcb_price(mbs_model.dl,mbs_model.th,mbs_model.ka,mbs_model.si,mbs_model.x0,mbs_model.T)
println("zcb price ",zcb_price)
println("Q at t=0 ",Q_f(0.0))

R_mc = get_R(mbs_model.t_sim,mbs_model.t_sim,mbs_model.T,mbs_model.r,mbs_model.neg_int_r_h)
R_f = ip.LinearInterpolation(mbs_model.t_sim,R_mc) # Q function
R_tm = R_f(tm) # Q value at tm
println("R at t=0 ",R_f(0.0))

#=

r0 = sum(mbs_model.x0)
if mbs_model.num_states == 1
    k0 = r0
else
    k0 = mbs_model.x0[end]
end
cases = [Dict(:a=>0.024, :b=>0.0, :gamma=>0.0, :T_asterisk=>1000.0, :k=>k0), # no prepayment
        Dict(:a=>0.024, :b=>0.75, :gamma=>0.0, :T_asterisk=>2.5, :k=>k0),
        Dict(:a=>0.024, :b=>0.0, :gamma=>10.0, :T_asterisk=>2.5, :k=>k0),
        Dict(:a=>0.024, :b=>0.75, :gamma=>10.0, :T_asterisk=>2.5, :k=>k0)] 

# m_curves = []
# for i in 1:length(cases)
#     update!(mbs_model; cases[i]...)
#     push!(m_curves,get_m_curve(mbs_model))
# end

# # Use PyPlot instead of Plots
# PyPlot.plot(mbs_model.t_sim[2:end],m_curves[1][2:end],label="No prep.")
# for i in 2:length(cases)
#     PyPlot.plot(mbs_model.t_sim[2:end],m_curves[i][2:end],label="b="*string(cases[i][:b])*" gamma="*string(cases[i][:gamma]))
# end
# PyPlot.title("r_0 = "*string(r0))
# PyPlot.legend()

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

#=
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
=#
=#