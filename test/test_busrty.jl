using DelaySSAToolkit
using Catalyst
using Statistics
using Test
using BenchmarkTools

reltol = 0.1
@parameters t
@variables X(t)
burst_sup = 30
a, b = [0.0282, 3.46]
rxs = [Reaction(a*b^i/(1+b)^(i+1),nothing,[X],nothing,[i]) for i in 1:burst_sup]
rxs = vcat(rxs)
@named rs = ReactionSystem(rxs,t,[X],[])

# convert the ReactionSystem to a JumpSystem
jumpsys = convert(JumpSystem, rs, combinatoric_ratelaws=false)

u0 = [0]
de_chan0 = [[]]
tf = 200.
tspan = (0,tf)
dprob = DiscreteProblem(jumpsys,u0,tspan)
τ = 130.

delay_trigger = Dict([Pair(i, [1=>fill(τ, i)]) for i in 1:burst_sup])
delay_complete = Dict(1=>[1=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)



algos = [DelayDirect(), DelayMNRM(), DelayRejection(), DelayDirectCR()]
timestamps = [10, 20, 50, 200]

bursty_mean(t) = a*b*min(t,τ)
bursty_var(t) = 2*a*b^2*min(t,τ) + a*b*min(t,τ)
samplesize = Int64(5e4)

function runbenchmark(algo)
    jprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(false, false))
    ensprob = EnsembleProblem(jprob)
    a=@benchmark solve(ensprob, SSAStepper(), EnsembleSerial(), trajectories=samplesize, saveat = timestamps)
    return a
end

function testalgo(algo_list,times)
    for algo in algo_list
        onealgotext=zeros(1,times)
        for i=1:times
            bench=runbenchmark(algo)
            onealgotext[i]=median(bench).time
        end
        push!(test_all_algo,onealgotext)
    end
end

times=2
algo_list = [DelayMNRM(), DelayRejection(), DelayDirectCR()]

test_all_algo=[]
@time testalgo(algo_list,times)

test_all_algo

#mean and std
#DelayMNRM()
mean(test_all_algo[1])
std(test_all_algo[1])
#DelayRejection()
mean(test_all_algo[2])
std(test_all_algo[2])
#DelayDirectCR()
mean(test_all_algo[3])
std(test_all_algo[3])

#=
jprob = DelayJumpProblem(jumpsys, dprob, algos[1], delayjumpset, de_chan0, save_positions=(false, false))
ensprob = EnsembleProblem(jprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleSerial(), trajectories=samplesize, saveat = timestamps)
a=@benchmark solve(ensprob, SSAStepper(), EnsembleSerial(), trajectories=samplesize, saveat = timestamps)
median(a).time
median(a)
median(a).gctime
=#