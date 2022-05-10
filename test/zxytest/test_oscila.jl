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
# alg = DelayDirectCR()
alg = DelayMNRM()
djprob = DelayJumpProblem(jumpsys, dprob, alg, delayjumpset, de_chan0)

sol = @benchmark solve(djprob, SSAStepper())

meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
function testalgo(algo_list)
    for algo in algo_list
        djprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0)
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
CSV.write("C:/Users/86158/Desktop/algotest/oscila_save_positions=default.csv",df)