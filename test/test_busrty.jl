using DelaySSAToolkit
using Catalyst
using Statistics
using Test
using BenchmarkTools
using DataFrames
using CSV
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

jprob = DelayJumpProblem(jumpsys, dprob, algos[1], delayjumpset, de_chan0, save_positions=(false, false))
a=@benchmark solve(jprob, SSAStepper(),saveat = timestamps)
mean(a)
median(a)
minimum(a)
maximum(a)
std(a)


meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
function testalgo(algo_list)
    for algo in algo_list
        jprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(false, false))
        a=@benchmark solve(jprob, SSAStepper(), saveat = timestamps)
        push!(meanlist,mean(a).time/1e9)
        push!(medianlist,median(a).time/1e9)
        push!(minlist,minimum(a).time/1e9)
        push!(maxlist,maximum(a).time/1e9)
        push!(stdlist,std(a).time/1e9)
    end
end

algo_list = [DelayMNRM(), DelayRejection(), DelayDirectCR()]

@time testalgo(algo_list)

meanlist
medianlist
minlist
maxlist
stdlist
df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("C:/Users/86158/Desktop/algotest/bursty_save_positions=f.csv",df)


using Printf
strmean = [@sprintf("%.5fs", yi) for yi in meanlist]
strstd = [@sprintf("%.3f", yi) for yi in stdlist]

sa=[string(round(mt,digits=7),"s") for mt in meanlist]
algo_name = ["MNRM", "Rejection", "DirectCR"]
using Plots
bar(algo_name,meanlist,legend=:false)
scatter!(algo_name, 0.1e-5 .+ meanlist , markeralpha=0, series_annotations=sa)