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

meanlist=[]
medianlist=[]
minlist=[]
maxlist=[]
stdlist=[]
for algo in algo_list
    djprob = DelayJumpProblem(jumpsys,dprob,algo,delayjumpset, de_chan0, save_positions=(false,false))
    a=@benchmark solve(djprob, SSAStepper())
    push!(meanlist,copy(mean(a).time/1e9))
    push!(medianlist,copy(median(a).time/1e9))
    push!(minlist,copy(minimum(a).time/1e9))
    push!(maxlist,copy(maximum(a).time/1e9))
    push!(stdlist,copy(std(a).time/1e9))
    print(djprob.aggregator," OK","\n")
end

df=DataFrame(mean=meanlist,median=medianlist,min=minlist,max=maxlist,std=stdlist)
CSV.write("test/zxytest/results/cascade_save_positions=f.csv",df)

df=CSV.read("test/zxytest/results/cascade_save_positions=f.csv",DataFrame)
medianlist=df.median

medianvalue=[string(round(mt,digits=6),"s") for mt in medianlist]

algo_name = ["DelayDirect","DelayRejection","DelayMNRM","DelayDirectCR"]
using Plots
p1=bar(algo_name,medianlist,legend=:false,title="cascade",ylabel="mediantime")
scatter!(algo_name, 0.0005 .+ medianlist , markeralpha=0, series_annotations=medianvalue)
savefig(p1,"test/zxytest/results/cascade_median.png")