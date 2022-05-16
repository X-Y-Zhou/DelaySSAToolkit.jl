using DelaySSAToolkit
using Test, DiffEqJump, StaticArrays
#σ_off: Gon -> Goff
#σ_on: Goff -> Gon
#ρ_on: Gon -> Gon + N, triggers N => τ 0
#d: P -> 0
# 1. Gon, 2. Goff, 3. P
params = [0.1, 0.1, 10., 0.1, 10., 50.]
σ_off, σ_on, ρ_on, d, τ, tf = params
rates = [σ_off, σ_on, ρ_on,  d]
react_stoich = [[1=>1],[2=>1],[1=>1],[3=>1]] 
net_stoich = [[1=>-1,2=>1],[1=>1,2=>-1],[3=>1],[3=>-1]]
mass_action_jump = MassActionJump(rates, react_stoich, net_stoich; scale_rates = false)
jumpset = JumpSet((),(),nothing,mass_action_jump)
delay_trigger = Dict(3=>[1=>τ])
delay_complete = Dict(1=>[3=>-1])
delay_affect! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1],i)
end
delay_interrupt = Dict(4=>delay_affect!)

delayjumpset = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
u0 =@SVector [0,1,0]
de_chan0 = [[]]
tspan = (0.,tf)
dprob = DiscreteProblem(u0, tspan)

djprob = DelayJumpProblem(dprob, DelayDirect(), jumpset, delayjumpset, de_chan0, save_positions = (true,true),save_delay_channel = true)
sol = @benchmark solve(djprob, SSAStepper())

algo_list =[DelayDirect(),DelayRejection(),DelayMNRM(),DelayDirectCR()]
meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
algo_list[1]

djprob = DelayJumpProblem(dprob, DelayMNRM(), jumpset, delayjumpset, de_chan0, save_positions = (true,true),save_delay_channel = true)
    a=@benchmark solve(djprob, SSAStepper())
    push!(meanlist,copy(mean(a).time/1e9))
    push!(medianlist,copy(median(a).time/1e9))
    push!(minlist,copy(minimum(a).time/1e9))
    push!(maxlist,copy(maximum(a).time/1e9))
    push!(stdlist,copy(std(a).time/1e9))
    print(djprob.aggregator," OK","\n")

1:length(integrator.de_chan[1])
1:1
length([])
@time testalgo(algo_list)

meanlist
medianlist
minlist
maxlist
stdlist
df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("C:/Users/86158/Desktop/algotest/save_delay_channel_save_positions=T.csv",df)