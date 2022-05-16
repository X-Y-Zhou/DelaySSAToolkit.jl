using DelaySSAToolkit, Catalyst
using Test
rn = @reaction_network begin
    1/(1+Y^2), 0 --> X
    1/(1+Y),   Y --> 0
end
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
# states(rn)
u0 = [0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
τ = 20.
dprob = DiscreteProblem(jumpsys, u0, tspan)
# jumpsys.dep_graph

delay_trigger_affect! = function (integrator, rng)
  append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 1=>-1], 2=>[2=>1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
algo_list =[DelayDirect(),DelayRejection(),DelayMNRM(),DelayDirectCR()]
for algo in algo_list
    djprob = DelayJumpProblem(jumpsys, dprob,algo,delayjumpset, de_chan0, save_positions=(false,false))
    a=@benchmark solve(djprob, SSAStepper())
    push!(meanlist,mean(a).time/1e9)
    push!(medianlist,median(a).time/1e9)
    push!(minlist,minimum(a).time/1e9)
    push!(maxlist,maximum(a).time/1e9)
    push!(stdlist,std(a).time/1e9)
    print(djprob.aggregator," OK","\n")
end

df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("test/zxytest/results/oscillation_save_positions=f.csv",df)

df=CSV.read("test/zxytest/results/oscillation_save_positions=f.csv",DataFrame)
medianlist=df.median

medianvalue=[string(round(mt,digits=6),"s") for mt in medianlist]

algo_name = ["DelayDirect","DelayRejection","DelayMNRM","DelayDirectCR"]
using Plots
p1=bar(algo_name,medianlist,legend=:false,title="oscillation",ylabel="mediantime")
scatter!(algo_name, 0.0003 .+ medianlist , markeralpha=0, series_annotations=medianvalue)
savefig(p1,"test/zxytest/results/oscillation_median.png")