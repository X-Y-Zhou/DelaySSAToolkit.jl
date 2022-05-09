using DiffEqJump, Catalyst, DelaySSAToolkit
rn = @reaction_network begin
   C, 0 --> Xₐ
   γ, Xₐ --> 0
   β, Xₐ --> Xᵢ
   γ, Xᵢ --> 0
end C γ β
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)

u0 = [0, 0]
tf = 30.
saveat = .1
C, γ, β = [2., 0.1, 0.5]
p = [C, γ, β]
tspan = (0.,tf)
dprob = DiscreteProblem(u0, tspan, p)

τ = 15.
delay_trigger_affect! = function (integrator, rng)
   append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>-1]) 
delay_affect! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1],i)
end
delay_interrupt = Dict(4=>delay_affect!) 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)

de_chan0 = [[]]
djprob = DelayJumpProblem(jumpsys, dprob, DelayMNRM(),  delaysets, de_chan0, save_positions=(false,false))

sol =@benchmark solve(djprob, SSAStepper(), seed=2, saveat =.1)

meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
function testalgo(algo_list)
    for algo in algo_list
        djprob = DelayJumpProblem(jumpsys, dprob,algo,delaysets, de_chan0, save_positions=(false,false))
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
CSV.write("C:/Users/86158/Desktop/algotest/dep_gr_save_positions=F.csv",df)