### Compares to the manually calcualted function ###
identical_networks_1 = Vector{Pair}()

rate1(u,p,t) = p[1]
rate2(u,p,t) = p[2]*u[1]
rate3(u,p,t) = p[3]*u[2]
rate4(u,p,t) = p[4]*u[2]
rate5(u,p,t) = p[5]*u[3]
rate6(u,p,t) = p[6]*u[3]
rate7(u,p,t) = p[7]*u[4]
rate8(u,p,t) = p[8]*u[4]
affect1!(integrator) = (integrator.u[1] += 1;)
affect2!(integrator) = (integrator.u[1] -= 1; integrator.u[2] += 1;)
affect3!(integrator) = (integrator.u[2] -= 1; integrator.u[1] += 1)
affect4!(integrator) = (integrator.u[2] -= 1; integrator.u[3] += 1)
affect5!(integrator) = (integrator.u[3] -= 1; integrator.u[2] += 1)
affect6!(integrator) = (integrator.u[3] -= 1; integrator.u[4] += 1)
affect7!(integrator) = (integrator.u[4] -= 1; integrator.u[3] += 1)
affect8!(integrator) = (integrator.u[4] -= 1;)
jump1 = ConstantRateJump(rate1,affect1!)
jump2 = ConstantRateJump(rate2,affect2!)
jump3 = ConstantRateJump(rate3,affect3!)
jump4 = ConstantRateJump(rate4,affect4!)
jump5 = ConstantRateJump(rate5,affect5!)
jump6 = ConstantRateJump(rate6,affect6!)
jump7 = ConstantRateJump(rate7,affect7!)
jump8 = ConstantRateJump(rate8,affect8!)
jumps = (jump1,jump2,jump3,jump4,jump5,jump6,jump7,jump8)
push!(identical_networks_1, reaction_networks_standard[6] => jumps)

rate1(u,p,t) = p[1]/10 + u[1]^p[3]/(u[1]^p[3] + p[2]^p[3])
rate2(u,p,t) = p[1]/10 + u[1]^p[3]/(u[1]^p[3] + p[2]^p[3])
rate3(u,p,t) = p[4]*u[1]*u[2]
rate4(u,p,t) = p[5]*u[3]
rate5(u,p,t) = p[6]*u[3]
rate6(u,p,t) = p[7]
rate7(u,p,t) = p[7]
rate8(u,p,t) = p[7]
affect1!(integrator) = (integrator.u[1] += 1;)
affect2!(integrator) = (integrator.u[2] += 1;)
affect3!(integrator) = (integrator.u[1] -= 1; integrator.u[2] -= 1; integrator.u[3] += 1)
affect4!(integrator) = (integrator.u[1] += 1; integrator.u[2] += 1; integrator.u[3] -= 1)
affect5!(integrator) = (integrator.u[3] -= 1; integrator.u[1] += 1;)
affect6!(integrator) = (integrator.u[1] -= 1;)
affect7!(integrator) = (integrator.u[2] -= 1;)
affect8!(integrator) = (integrator.u[3] -= 1;)
jump1 = ConstantRateJump(rate1,affect1!)
jump2 = ConstantRateJump(rate2,affect2!)
jump3 = ConstantRateJump(rate3,affect3!)
jump4 = ConstantRateJump(rate4,affect4!)
jump5 = ConstantRateJump(rate5,affect5!)
jump6 = ConstantRateJump(rate6,affect6!)
jump7 = ConstantRateJump(rate7,affect7!)
jump8 = ConstantRateJump(rate8,affect8!)
jumps = (jump1,jump2,jump3,jump4,jump5,jump6,jump7,jump8)
push!(identical_networks_1, reaction_networks_hill[7] => jumps)

rate1(u,p,t) = p[1]*binomial(u[1],1)
rate2(u,p,t) = p[2]*binomial(u[2],2)
rate3(u,p,t) = p[3]*binomial(u[2],2)
rate4(u,p,t) = p[4]*binomial(u[3],3)
rate5(u,p,t) = p[5]*binomial(u[3],3)
rate6(u,p,t) = p[6]*binomial(u[4],4)
affect1!(integrator) = (integrator.u[1] -=1 ; integrator.u[2] +=2 ;)
affect2!(integrator) = (integrator.u[2] -=2 ; integrator.u[1] +=1 ;)
affect3!(integrator) = (integrator.u[2] -=2 ; integrator.u[3] +=3 ;)
affect4!(integrator) = (integrator.u[3] -=3 ; integrator.u[2] +=2 ;)
affect5!(integrator) = (integrator.u[3] -=3 ; integrator.u[4] +=4 ;)
affect6!(integrator) = (integrator.u[4] -=4 ; integrator.u[3] +=3 ;)
jump1 = ConstantRateJump(rate1,affect1!)
jump2 = ConstantRateJump(rate2,affect2!)
jump3 = ConstantRateJump(rate3,affect3!)
jump4 = ConstantRateJump(rate4,affect4!)
jump5 = ConstantRateJump(rate5,affect5!)
jump6 = ConstantRateJump(rate6,affect6!)
jumps = (jump1,jump2,jump3,jump4,jump5,jump6)
push!(identical_networks_1, reaction_networks_constraint[5] => jumps)

for networks in identical_networks_1
    println("HERE")
    for factor in [1e-2, 1e-1, 1e0, 1e1], repeat = 1:5
        u0 = rand(1:Int64(factor*100),length(networks[1].states))
        p = factor*rand(length(networks[1].ps))
        prob1 = JumpProblem(networks[1],DiscreteProblem(networks[1],u0,(0.,10000.),p),Direct())
        sol1 = solve(prob1,SSAStepper())
        prob2 = JumpProblem(DiscreteProblem(u0,(0.,10000.),p),Direct(),networks[2]...)
        sol2 = solve(prob2,SSAStepper())
        for i = 1:length(u0)
            vals1 = getindex.(sol1.u,i);
            vals2 = getindex.(sol1.u,i);
            @test 0.8 < mean(vals1)/mean(vals2) < 1.25
            @test 0.8 < std(vals1)/std(vals2) < 1.25
        end
    end
end
