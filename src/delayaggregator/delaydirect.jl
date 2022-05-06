mutable struct DelayDirectJumpAggregation{T,S,F1,F2,RNG,IType} <: AbstractDSSAJumpAggregator
    next_jump::Int
    prev_jump::Int
    next_jump_time::T
    end_time::T
    cur_rates::Vector{T} # cur_rates here is the cumsum of rates
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
    rng::RNG
    next_delay::Union{Nothing,Vector{Int}}
    num_next_delay::Union{Nothing,Vector{Int}}
    time_to_next_jump::T
    shadow_integrator::IType
end
mutable struct Shadow_Integrator{uType,chanType}
    u::uType
    de_chan::chanType
    delayjumpsets::DelayJumpSet
end
function DelayDirectJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; u0, kwargs...) where {T,S,F1,F2,RNG}
    ttnj = zero(et)
    shadow_integrator = Shadow_Integrator{typeof(u0),Vector{Vector{T}}}(copy(u0), [Vector{T}()], DelayJumpSet(Dict(), Dict(), Dict()))
    nd = nothing
    nnd = [1] # in the Direct Method the number of next delay equals always 1
    DelayDirectJumpAggregation{T,S,F1,F2,RNG,typeof(shadow_integrator)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, nd, nnd, ttnj, shadow_integrator)
end

function aggregate(aggregator::DelayDirect, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(DelayDirectJumpAggregation, u, p, t, end_time, ma_jumps,
        rates, affects!, save_positions, rng; u0=copy(u), kwargs...)
end

function initialize!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    p.shadow_integrator.delayjumpsets = integrator.delayjumpsets
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

@inline function execute_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    update_state_delay_Direct!(p, integrator, integrator.u, t)
    nothing
end

function generate_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    generate_delta!(p, integrator, params, t)
    @fastmath p.next_jump_time = t + p.time_to_next_jump
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # 
    nothing
end
"""
    Create delta based on the shawdow variable u_shadow
"""
@inline function generate_delta!(p::DelayDirectJumpAggregation, integrator, params, t)
    @unpack delay_complete = integrator.delayjumpsets
    p.shadow_integrator.u = copy(integrator.u) #TODO
    p.shadow_integrator.de_chan = deepcopy(integrator.de_chan) #TODO

    direct_algo!(p, integrator, t)
end

function direct_algo!(p, integrator, t)
    fill_cum_rates_and_sum!(p, p.shadow_integrator.u, params, t)
    r1 = rand(p.rng)
    if isempty(reduce(vcat, integrator.de_chan))
        ttnj = -log(r1) / p.sum_rate
        ttnj_last = ttnj
    else
        T1, T2 = create_Tstruct(integrator.de_chan) # TODO change the order, maybe create T1[2] first check if F < r1, then create Tstruct1
        prepend!(T1, zero(t))
        append!(T1, Inf)
        i = 1
        aₜ = p.sum_rate * T1[2]
        F = one(t) - exp(-aₜ)
        aₜ_ = zero(aₜ)
        sum_rate_ = calculate_sum_rate(p, p.shadow_integrator.u, params, t)
        while F < r1
            p.next_delay = [T2[i]]

            # TODO change the de_chan step by step
            shift_delay_channel!(p.shadow_integrator.de_chan, T1[i+1] - T1[i])
            update_delay_channel!(p.shadow_integrator.de_chan)
            update_delay_complete!(p, p.shadow_integrator)

            sum_rate_ = calculate_sum_rate(p, p.shadow_integrator.u, params, t + T1[i+1])

            aₜ_ = copy(aₜ) # backup aₜ
            aₜ += sum_rate_ * (T1[i+2] - T1[i+1])
            F = one(t) - exp(-aₜ)
            i += 1
        end
        p.sum_rate = sum_rate_
        ttnj_last = (-log(one(t) - r1) - aₜ_) / p.sum_rate
        # TODO change
        ttnj = T1[i] + ttnj_last
    end
    T1_last, T2_last = create_Tstruct(p.shadow_integrator.de_chan)

    shift_delay_channel!(p.shadow_integrator.de_chan, ttnj_last)
    update_delay_channel!(p.shadow_integrator.de_chan)

    # in case the last ttnj also change the state
    update_state_final_jump!(p, p.shadow_integrator, ttnj_last, T1_last, T2_last)

    fill_cum_rates_and_sum!(p, p.shadow_integrator.u, params, t + ttnj)
    p.time_to_next_jump = ttnj
end

@inbounds function update_state_final_jump!(p, integrator, tgap, T1, T2)
    idx = count(x -> x <= tgap, T1)
    if !isempty(idx)
        for i in 1:idx
            p.next_delay = [T2[i]]
            update_delay_complete!(p, integrator)
        end
    end
    nothing
end

function update_delay_chan_state_at_tstop!(p, integrator, tgap, T1, T2)
    idx = count(x -> x <= tgap, T1)
    # T1copy = copy(T1)
    # prepend!(T1copy, zero(tgap))
    # diff_T1 = diff(T1copy)
    if idx > 0
        ttnj_last = tgap - T1[idx]
        prev_T1 = 0
        for i in 1:idx
            p.next_delay = [T2[i]]
            shift_delay_channel!(integrator.de_chan, T1[i]- prev_T1)
            update_delay_channel!(integrator.de_chan)
            update_delay_complete!(p, integrator)
            prev_T1 = T1[i]
        end
    else
        ttnj_last = tgap
    end
    deleteat!(T1,1:idx)
    deleteat!(T2,1:idx)
    
    shift_delay_channel!(integrator.de_chan, ttnj_last)
    update_delay_channel!(integrator.de_chan)
    nothing
end

"""
    function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)

This function modifies `p.cur_rates` and `p.sum_rate`
"""
function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = p.cur_rates
    # mass action rates
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in 1:length(p.rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end

    @inbounds sum_rate = cur_rates[end]
    p.sum_rate = sum_rate
    nothing
end

function calculate_sum_rate(p::DelayDirectJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = zeros(length(p.cur_rates))
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in 1:length(p.rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end
    cur_rates[end]
end
"""
    function create_Tstruct(de_chan::Vector{Vector{T}})

calculate `Tstruct` according to the de_chan.
"""
function create_Tstruct(de_chan::Vector{Vector{T}}) where {T}
    N = sum(length.(de_chan))
    Tstruct1 = Vector{T}(undef, N)
    Tstruct2 = Vector{Int64}(undef, N)
    k = 1
    @inbounds while k <= N
        for i in eachindex(de_chan)
            for j in eachindex(de_chan[i])
                Tstruct1[k] = de_chan[i][j]
                Tstruct2[k] = i
                k += 1
            end
        end
    end
    vecorder = sortperm(Tstruct1)
    Tstruct1 = Tstruct1[vecorder]
    Tstruct2 = Tstruct2[vecorder]
    Tstruct1, Tstruct2
end

@inline function update_state_delay_Direct!(p::DelayDirectJumpAggregation, integrator, u, t)
    @unpack ma_jumps, next_jump, time_to_next_jump = p
    @unpack delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets

    integrator.u = copy(p.shadow_integrator.u)
    integrator.de_chan = deepcopy(p.shadow_integrator.de_chan) #TODO 

    num_ma_rates = get_num_majumps(ma_jumps)
    if next_jump <= num_ma_rates # if the next jump is a mass action jump
        if u isa SVector
            integrator.u = executerx(u, next_jump, ma_jumps)
        else
            @inbounds executerx!(integrator.u, next_jump, ma_jumps)
        end
    else
        idx = next_jump - num_ma_rates
        @inbounds p.affects![idx](integrator)
    end

    if next_jump in delay_interrupt_set
        update_delay_interrupt!(p, integrator)
    elseif next_jump in delay_trigger_set
        update_delay_trigger!(p, integrator)
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end

