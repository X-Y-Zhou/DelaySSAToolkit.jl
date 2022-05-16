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

algo_list = [DelayDirect(),DelayRejection(),DelayMNRM(),DelayDirectCR()]

jprob = DelayJumpProblem(jumpsys, dprob, algo_list[1], delayjumpset, de_chan0, save_positions=(false, false))
a=@benchmark solve(jprob, SSAStepper())

jprob = DelayJumpProblem(jumpsys, dprob, algo_list[2], delayjumpset, de_chan0, save_positions=(false, false))
a=@benchmark solve(jprob, SSAStepper())

jprob = DelayJumpProblem(jumpsys, dprob, algo_list[3], delayjumpset, de_chan0, save_positions=(false, false))
a=@benchmark solve(jprob, SSAStepper())

jprob = DelayJumpProblem(jumpsys, dprob,algo_list[4] , delayjumpset, de_chan0, save_positions=(false, false))
a=@benchmark solve(jprob, SSAStepper())


meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
function testalgo(algo_list)
    for algo in algo_list
        #de_chan0 = [[]]
        jprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(false, false))
        a=@benchmark solve(jprob, SSAStepper())
        push!(meanlist,mean(a).time/1e9)
        push!(medianlist,median(a).time/1e9)
        push!(minlist,minimum(a).time/1e9)
        push!(maxlist,maximum(a).time/1e9)
        push!(stdlist,std(a).time/1e9)
    end
end

@time testalgo(algo_list)

df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("test/zxytest/results/bursty_save_positions=f.csv",df)

df=CSV.read("C:/Users/86158/Desktop/algotest/bursty_save_positions=f.csv",DataFrame)

medianlist=df.median

algo_name = ["Direct","Rejection","MNRM","DirectCR"]

meanvalue=[string(round(mt,digits=8),"s") for mt in meanlist]
medianvalue=[string(round(mt,digits=8),"s") for mt in medianlist]

using Plots
bar(algo_name,meanlist,legend=:false,title="bursty",ylabel="meantime")
scatter!(algo_name, 0.3e-5 .+ meanlist , markeralpha=0, series_annotations=meanvalue)

bar(algo_name,medianlist,legend=:false,title="bursty",ylabel="mediantime")
scatter!(algo_name, 0.3e-5 .+ medianlist , markeralpha=0, series_annotations=medianvalue)

