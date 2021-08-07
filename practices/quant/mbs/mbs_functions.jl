
import Distributions; distr = Distributions
import Statistics; stats = Statistics
import NumericalIntegration; ni = NumericalIntegration
import Random



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


# integrate r from u to T where 0<=u<T
"""
`l`: shift function
"""
function get_int_r(t_sim,x_u,l=nothing)
    num_steps = length(t_sim)
    num_states = length(x_u)
    num_paths = size(x_u[1],1)
    int_r = zeros(num_paths, num_steps)
    for sv_idx in 1:num_states
        for jdx in 1:num_steps-1 # exclude the last
            for idx in 1:num_paths
                int_r[idx,jdx] += ni.integrate(@view(t_sim[jdx:end]), @view(x_u[sv_idx][idx,jdx:end]))
            end
        end
    end
    return int_r    
end

function get_r_at_t(t_sim,t,x_u,l=nothing)
    num_states = length(x_u)
    num_paths = size(x_u[1],1)
    t_ind = findfirst(x->x==t, t_sim)
    r = zeros(num_paths)
    for sv_idx in 1:num_states
        r += x_u[sv_idx][:,t_ind]
    end
    return r
end

function get_r(x_u,l=nothing)
    num_states = length(x_u)
    num_paths = size(x_u[1],1)
    num_steps = size(x_u[1],2)
    r = zeros(size(x_u[1]))
    for sv_idx in 1:num_states
        for jdx in 1:num_steps
            for idx in 1:num_paths
                r[idx,jdx] += x_u[sv_idx][idx,jdx]
            end
        end
    end
    return r
end


function get_inter_h0_1!(inter_h0_1,t_sim,T_asterisk)
    h0(t) = min(t,T_asterisk)
    for ind in 1:length(t_sim)-1
        if ind == length(t_sim) # skip the last
            continue
        end
        h0_array = h0.(@view t_sim[ind:end])
        inter_h0_1[ind] = ni.integrate(@view(t_sim[ind:end]), h0_array)
    end
    return inter_h0_1
end

function get_inter_h_r!(inter_h0_2,t_sim,k,r)
    @. inter_h0_2 = max(0.0,k -r) # k-r(s)
    num_paths = size(inter_h0_2,1)
    # compute ∫(k-r(s))ds 
    for jdx in 1:length(t_sim)-1 # exclude the last
        for idx in 1:num_paths
            inter_h0_2[idx,jdx] += ni.integrate(@view(t_sim[jdx:end]), @view(inter_h0_2[idx,jdx:end]))
        end
    end
    return inter_h0_2
end

function get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
    get_inter_h0_1!(inter_h0_1,t_sim,T_asterisk)
    get_inter_h_r!(inter_h0_2,t_sim,k,r)
    return inter_h0_1, inter_h0_2
end

function get_inter_h0(t_sim,k,r,T_asterisk)
    num_steps = size(r,2) - 1
    inter_h0_1 = zeros(num_steps+1)
    inter_h0_2 = zeros(size(r))
    return get_inter_h0!(inter_h0_1, inter_h0_2,t_sim,k,r,T_asterisk)
end

function get_int_h!(int_h, ab, gamma, inter_h0_1, inter_h0_2)
    int_h .= gamma * inter_h0_2
    tmp = ab*inter_h0_1
    for idx in 1:size(int_h,1) # for each path
        int_h[idx,:] .+= tmp
    end
end

"""
Compute E[rᵤe^(∫ₜᵘ-(rₛ+hₛ)ds)]
"""
function get_R(t_sim,t,u,r,neg_int_r_h)
    t_ind = findfirst(x->x==t,t_sim)
    u_ind = findfirst(x->x==u,t_sim)
    r_u = r[:,u_ind]
    neg_int_t_u = neg_int_r_h[:,t_ind] - neg_int_r_h[:,u_ind] # ∫ₜᵘ-(rₛ+hₛ)ds
    exp_neg_int = exp.(neg_int_t_u)
    res = stats.mean(r_u.*exp_neg_int)
    return res
end

"""
Compute E[e^(∫ₜᵘ-(rₛ+hₛ)ds)]
"""
function get_Q(t_sim,t,u,neg_int_r_h)
    t_ind = findfirst(x->x==t,t_sim)
    u_ind = findfirst(x->x==u,t_sim)
    neg_int_t_u = neg_int_r_h[:,t_ind] - neg_int_r_h[:,u_ind] # ∫ₜᵘ-(rₛ+hₛ)ds
    exp_neg_int = exp.(neg_int_t_u)
    res = stats.mean(exp_neg_int)
    return res
end

function m_equation(m0,t_sim,T,r,neg_int_r_h)
    t_space = @view(t_sim[0.0.<=t_sim.<=T])
    num_integrand = zeros(size(t_space))
    den_integrand = zeros(size(t_space))
    for (ind,u) in enumerate(t_space)
        if u == T
            continue # skip the last not to comput R and Q
        end
        R_0_u = get_R(t_space,0.0,u,r,neg_int_r_h)
        Q_0_u = get_Q(t_space,0.0,u,neg_int_r_h)
        factor = 1.0 - exp(-m0*(T-u))
        num_integrand[ind] = factor*R_0_u
        den_integrand[ind] = factor*Q_0_u
    end
    # @show num_integrand
    # @show den_integrand
    res = ni.integrate(t_space, num_integrand)
    res /= ni.integrate(t_space, den_integrand)
    return res
end

# fixed-iteration scheme
function get_m_fi(m0,t_sim,T,r,neg_int_r_h;max_iter=100,abs_tol=1e-6)
    mi = m0
    mp = m0
    for i in 1:max_iter
        # @show i # to check how many steps are necessary
        mi = m_equation(mi,t_sim,T,r,neg_int_r_h)
        if abs(mi-mp) <= abs_tol
            break
        end
        mp = mi
    end
    return mi
end

# simulations
function simulate_x(num_paths,num_steps,ka,si,th,s,T,x0)
    Random.seed!(2) # make sure the same seeding is used
    x_u = [] # x value at each step
    for i in 1:length(x0)
        push!(x_u, get_CIR_sample_steps(num_paths,num_steps,ka[i],si[i],th[i],s,T,x0[i]))
    end
    return x_u
end