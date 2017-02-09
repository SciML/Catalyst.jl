# Not actually used: just for comparison

using Gillespie
using Gadfly

function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end

x0 = [999,1,0]
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/1000.0,0.01]
tf = 250.0
srand(1234)

nums = Int[]
@time for i in 1:1000
  result = ssa(x0,F,nu,parms,tf)
  data = ssa_data(result)
  push!(nums,data[:x3][end])
end
mean(nums)

result = ssa(x0,F,nu,parms,tf)
data = ssa_data(result)

p=plot(data,
  layer(x="time",y="x1",Geom.step,Theme(default_color=color("red"))),
  layer(x="time",y="x2",Geom.step,Theme(default_color=color("blue"))),
  layer(x="time",y="x3",Geom.step,Theme(default_color=color("green"))),
  Guide.xlabel("Time"),
  Guide.ylabel("Number"),
  Guide.manual_color_key("Population",
                            ["S", "I", "R"],
                            ["red", "blue", "green"]),
  Guide.title("SIR epidemiological model"))
