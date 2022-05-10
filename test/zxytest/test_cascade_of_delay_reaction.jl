using Catalyst, DelaySSAToolkit
# 0 -->A1, => A2, ..., => A5 => 0 with delay 1 for each cascade delay

rn = @reaction_network begin
    5, 0 --> A1
end

delay_trigger_affect! = []
chain_len = 10
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], 1.)
end


u0 = zeros(chain_len)
de_chan0 = []
for _ in 1:chain_len
    push!(de_chan0, [])
end
tspan = (0.,50.)
delay_complete_affect! = []
for i in 1:chain_len-1
    push!(delay_complete_affect!, function (integrator, rng)
        integrator.u[i] -= 1 # A_prev minus 1
        integrator.u[i+1] += 1 # A plus 1
        append!(integrator.de_chan[i+1], 1.) # add to the delay channel
    end
    )
end
push!(delay_complete_affect!, function (integrator, rng)
    integrator.u[chain_len] -= 1 # A_prev minus 1
end)

delay_trigger = Dict(Pair(1, delay_trigger_affect!))
delay_complete = Dict(i=>delay_complete_affect![i] for i in 1:chain_len)
delay_interrupt = Dict()


delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
dprob = DiscreteProblem(jumpsys, u0, tspan)

algo_list =[DelayDirect(),DelayRejection(),DelayMNRM(),DelayDirectCR()]

djprob = DelayJumpProblem(jumpsys,dprob,algo_list[1],delayjumpset, de_chan0, save_positions=(false,false))
a=@benchmark solve(djprob, SSAStepper(), saveat =.1)

djprob = DelayJumpProblem(jumpsys,dprob,algo_list[2],delayjumpset, de_chan0, save_positions=(false,false))
a=@benchmark solve(djprob, SSAStepper(), saveat =.1)

djprob = DelayJumpProblem(jumpsys,dprob,algo_list[3],delayjumpset, de_chan0, save_positions=(false,false))
a=@benchmark solve(djprob, SSAStepper(), saveat =.1)

djprob = DelayJumpProblem(jumpsys,dprob,algo_list[4],delayjumpset, de_chan0, save_positions=(false,false))
a=@benchmark solve(djprob, SSAStepper(), saveat =.1)

function testalgo(algo_list)
    for algo in algo_list
        djprob = DelayJumpProblem(jumpsys,dprob,algo,delayjumpset, de_chan0, save_positions=(false,false))
        a=@benchmark solve(djprob, SSAStepper(), saveat =.1)
        push!(meanlist,copy(mean(a).time/1e9))
        push!(medianlist,copy(median(a).time/1e9))
        push!(minlist,copy(minimum(a).time/1e9))
        push!(maxlist,copy(maximum(a).time/1e9))
        push!(stdlist,copy(std(a).time/1e9))
    end
end
meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
@time testalgo(algo_list)

df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("C:/Users/86158/Desktop/algotest/cascade_save_positions=F.csv",df)

sa=[string(round(mt,digits=6),"s") for mt in meanlist]
algo_name = ["Direct","Rejection","MNRM","DirectCR"]
using Plots
bar(algo_name,meanlist,legend=:false,title="cascade",ylabel="meantime")
scatter!(algo_name, 0.001 .+ meanlist , markeralpha=0, series_annotations=sa)
