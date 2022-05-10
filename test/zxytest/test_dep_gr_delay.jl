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

algo_list =[DelayDirect(),DelayRejection(),DelayMNRM(),DelayDirectCR()]

djprob = DelayJumpProblem(jumpsys, dprob, algo_list[1],  delaysets, de_chan0, save_positions=(false,false))
sol =@benchmark solve(djprob, SSAStepper(),saveat =.1)

djprob = DelayJumpProblem(jumpsys, dprob,algo_list[2],  delaysets, de_chan0, save_positions=(false,false))
sol =@benchmark solve(djprob, SSAStepper(),saveat =.1)

djprob = DelayJumpProblem(jumpsys, dprob,algo_list[3],  delaysets, de_chan0, save_positions=(false,false))
sol =@benchmark solve(djprob, SSAStepper(),saveat =.1)

djprob = DelayJumpProblem(jumpsys, dprob,algo_list[4],  delaysets, de_chan0, save_positions=(false,false))
sol =@benchmark solve(djprob, SSAStepper(),saveat =.1)

meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
function testalgo(algo_list)
    for algo in algo_list
        djprob = DelayJumpProblem(jumpsys, dprob,algo,delaysets, de_chan0, save_positions=(false,false))
        a=@benchmark solve(djprob, SSAStepper(), saveat =.1)
        push!(meanlist,mean(a).time/1e9)
        push!(medianlist,median(a).time/1e9)
        push!(minlist,minimum(a).time/1e9)
        push!(maxlist,maximum(a).time/1e9)
        push!(stdlist,std(a).time/1e9)
    end
end

@time testalgo(algo_list)

df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("C:/Users/86158/Desktop/algotest/dep_gr_save_positions=F.csv",df)

df=CSV.read("C:/Users/86158/Desktop/algotest/dep_gr_save_positions=F.csv",DataFrame)
medianlist=df.median

sa=[string(round(mt,digits=6),"s") for mt in meanlist]
medianvalue=[string(round(mt,digits=6),"s") for mt in medianlist]

algo_name = ["Direct","Rejection","MNRM","DirectCR"]
using Plots
bar(algo_name,meanlist,legend=:false,title="delay_degradation",ylabel="meantime")
scatter!(algo_name, 0.00015 .+ meanlist , markeralpha=0, series_annotations=sa)

bar(algo_name,medianlist,legend=:false,title="delay_degradation",ylabel="mediantime")
scatter!(algo_name, 0.00015 .+ medianlist , markeralpha=0, series_annotations=medianvalue)