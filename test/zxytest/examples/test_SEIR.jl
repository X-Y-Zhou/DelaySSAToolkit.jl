using Catalyst
using Test
using DelaySSAToolkit
rn = @reaction_network begin
    ρ, S+I --> E+I
    1/I, I --> R
    r, I --> R
end ρ r

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)

u0 = [999,1,0,0] # S, I, E, R
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2] # parameters for ρ, r
τ = 20.
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)

delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_interrupt = Dict()
delay_complete = Dict(1=>[2=>1, 3=>-1])
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

de_chan0 = [[]]

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
CSV.write("test/zxytest/results/SEIR_save_positions=f.csv",df)

df=CSV.read("test/zxytest/results/SEIR_save_positions=f.csv",DataFrame)
medianlist=df.median

medianvalue=[string(round(mt,digits=6),"s") for mt in medianlist]

algo_name = ["DelayDirect","DelayRejection","DelayMNRM","DelayDirectCR"]
using Plots
p1=bar(algo_name,medianlist,legend=:false,title="SEIR",ylabel="mediantime")
scatter!(algo_name, 0.0002 .+ medianlist , markeralpha=0, series_annotations=medianvalue)
savefig(p1,"test/zxytest/results/SEIR_median.png")

#=
function runbenchmark(algo)
    djprob = DelayJumpProblem(dprob, algo, jumpset, delayjumpset, de_chan0, save_positions = (true,true),save_delay_channel = true)
    a=@benchmark solve(djprob, SSAStepper())
    #a=@bprofile solve(jprob, SSAStepper(), saveat = timestamps)
    return a
end

function testalgo(algo_list,times)
    for algo in algo_list
        onealgotext=zeros(1,times)
        for i=1:times
            bench=runbenchmark(algo)
            onealgotext[i]=median(bench).time/1e3
        end
        push!(test_all_algo,onealgotext)
    end
end


times=5
algo_list = [DelayMNRM(), DelayRejection(), DelayDirectCR()]

test_all_algo=[]
@time testalgo(algo_list,times)

test_all_algo
length(test_all_algo[1])
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

algo_name = ["MNRM", "Rejection", "DirectCR"]
meanlist=[mean(test_all_algo[i]) for i=1:length(test_all_algo)]
stdlist=[std(test_all_algo[i]) for i=1:length(test_all_algo)]

using Printf
strmean = [@sprintf("%.3fus", yi) for yi in meanlist]
strstd = [@sprintf("%.3f", yi) for yi in stdlist]

using Plots
p1=bar(algo_name,meanlist,ylims=(16,22),text=strmean,bar_width=0.1,
       label="mean of time",text_position="outside")
p2=bar(algo_name,stdlist,ylims=(0.15,0.5),text=strstd,bar_width=0.1,
       label="std of time")
plot(p1, p2, layout = (2,1))
=#