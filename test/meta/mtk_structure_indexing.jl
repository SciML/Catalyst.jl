# Fetch packages
using Catalyst, OrdinaryDiffEq, StochasticDiffEq, JumpProcesses, NonlinearSolve, Plots
using Test

# Create model, problems, and solutions.
begin
    model = complete(@reaction_network begin
        @observables X2 ~ 2X
        (kp,kd), 0 <--> X
        (k1,k2), X <--> Y
    end)
    @unpack X, Y, kp, kd, k1, k2 = model

    u0_vals = [X => 1, Y => 0]
    tspan = (0.0, 10.0)
    p_vals = [kp => 1.0, kd => 0.2, k1 => 1.0, k2 => 2.0]

    oprob = ODEProblem(model, u0_vals, tspan, p_vals)
    sprob = SDEProblem(model,u0_vals, tspan, p_vals)
    dprob = DiscreteProblem(model, u0_vals, tspan, p_vals)
    jprob = JumpProblem(model, dprob, Direct())
    nprob = NonlinearProblem(model, u0_vals, p_vals)
    problems = [oprob, sprob, dprob, jprob, nprob]
    
    osol = solve(oprob, Tsit5())
    ssol = solve(sprob, ImplicitEM())
    jsol = solve(jprob, SSAStepper())
    nsol = solve(nprob, NewtonRaphson())
    sols = [osol, ssol, jsol, nsol]
end

# Tests index updating.
let 
    for prob in deepcopy(problems)
        # Get u values.
        @test prob[X] == prob[model.X] == prob[:X] == 1
        @test prob[X2] == prob[model.X2] == prob[:X2] == 2
        @test prob[[X2,Y]] == prob[[model.X2,model.Y]] == prob[[:X2,:Y]] == [2, 2]
        @test prob[(X2,Y)] == prob[(model.X2,model.Y)] == prob[(:X2,:Y)] == (2, 2)
        @test getu(prob, X)(prob) == getu(prob, model.X)(prob) == getu(prob, :X)(prob) == 1   
        @test getu(prob, X2)(prob) == getu(prob, model.X2)(prob) == getu(prob, :X2)(prob) == 1      

        # Set u values.
        prob[X] = 2
        @test prob[X] == 2
        prob[model.X] = 3
        @test prob[X] == 3
        prob[:X] = 4
        @test prob[X] == 4
        setu(prob, X)(prob, 5)
        @test prob[X] == 5
        setu(prob, model.X)(prob, 6)
        @test prob[X] == 6
        setu(prob, :X)(prob, 7)
        @test prob[X] == 7

        # Get p values.
        @test getp(prob, kp)(prob) == getp(prob, model.kp)(prob) == getp(prob, :kp)(prob) == 1.0
        @test prob.ps[kp] == prob.ps[model.kp] == prob.ps[:kp] == 1.0    
        
        # Set p values.
        setp(prob, kp)(prob, 2.0)
        @test prob[kp] == 2.0
        setp(prob, model.kp)(prob, 3.0)
        @test prob[kp] == 3.0
        setp(prob, :kp)(prob, 4.0)
        @test prob[kp] == 4.0
        prob.ps[kp] = 5.0
        @test prob[kp] == .0
        prob.ps[model.kp] = 6.0
        @test prob[kp] == 6.0
        prob.ps[:kp] = 7.0
        @test prob[kp] == 7.0
    end
end

# Test remake function.
let 
    for prob in deepcopy(probs)
        # Remake for all u0s.
        @test remake(prob; u0 = [X => 2, Y => 3]).u0 == [2, 3]
        @test remake(prob; u0 = [model.X => 4, model.Y => 5]).u0 == [4, 5]
        @test remake(prob; u0 = [:X => 6, :Y => 7]).u0 == [6, 7]

        # Remake for some u0s.
        @test remake(prob; u0 = [Y => 8]).u0 == [1, 8]
        @test remake(prob; u0 = [model.Y => 9]).u0 == [1, 9]
        @test remake(prob; u0 = [:Y => 10]).u0 == [1, 10]

        # Remake for all ps.
        @test remake(prob; p = [kp => 1.0, kd => 2.0, k1 => 3.0, k2 => 4.0]).p == [1.0, 2.0, 3.0, 4.0]
        @test remake(prob; p = [model.kp => 5.0, model.kd => 6.0, model.k1 => 7.0, model.k2 => 8.0]).p == [5.0, 6.0, 7.0, 8.0]
        @test remake(prob; p = [:kp => 9.0, :kd => 10.0, :k1 => 11.0, :k2 => 12.0]).p  == [9.0, 10.0, 11.0, 12.0]

        # Remake for some ps.
        @test remake(prob; p = [k2 => 13.0]).p == [1.0, 0.2, 1.0, 13.0]
        @test remake(prob; p = [model.k2 => 14.0]).p == [1.0, 0.2, 1.0, 14.0]
        @test remake(prob; p = [:k2 => 15.0]).p  == [1.0, 0.2, 1.0, 15.0]
    end
end

# Test integrator indexing.
let 
    for integrator in init.(deepcopy(problems))
        # Get u values.
        @test integrator[X] == integrator[model.X] == integrator[:X] == 1
        @test integrator[X2] == integrator[model.X2] == integrator[:X2] == 2
        @test integrator[[X,Y]] == integrator[[model.X,model.Y]] == integrator[[:X,:Y]] == [1, 2]
        @test integrator[[X2,Y]] == integrator[[model.X2,model.Y]] == integrator[[:X2,:Y]] == [2, 2]
        @test integrator[(X,Y)] == integrator[(model.X,model.Y)] == integrator[(:X,:Y)] == (1, 2)   
        @test integrator[(X2,Y)] == integrator[(model.X2,model.Y)] == integrator[(:X2,:Y)] == (2, 2)   
        @test getu(integrator, X)(integrator) == getu(integrator, model.X)(integrator) == getu(integrator, :X)(integrator) == 1       
        @test getu(integrator, X2)(integrator) == getu(integrator, model.X2)(integrator) == getu(integrator, :X2)(integrator) == 2         

        # Set u values.
        integrator[X] = 2
        @test integrator[X] == 2
        integrator[model.X] = 3
        @test integrator[X] == 3
        integrator[:X] = 4
        @test integrator[X] == 4
        setu(integrator, X)(integrator, 5)
        @test integrator[X] == 5
        setu(integrator, model.X)(integrator, 6)
        @test integrator[X] == 6
        setu(integrator, :X)(integrator, 7)
        @test integrator[X] == 7

        # Get p values.
        @test getp(integrator, kp)(integrator) == getp(integrator, model.kp)(integrator) == getp(integrator, :kp)(integrator) == 1.0
        @test integrator.ps[kp] == integrator.ps[model.kp] == integrator.ps[:kp] == 1.0    

        # Set p values.
        setp(integrator, kp)(integrator, 2.0)
        @test integrator[kp] == 2.0
        setp(integrator, model.kp)(integrator, 3.0)
        @test integrator[kp] == 3.0
        setp(integrator, :kp)(integrator, 4.0)integrator
        @test integrator[kp] == 4.0
        integrator.ps[kp] = 5.0
        @test integrator[kp] == .0
        integrator.ps[model.kp] = 6.0
        @test integrator[kp] == 6.0
        integrator.ps[:kp] = 7.0
        @test integrator[kp] == 7.0
    end
end

# Test solve's save_idxs argument.
let 
    @test length(solve(oprob, Tsit5(); save_idxs=[X]).u[1]) == 1
    @test length(solve(sprob, ImplicitEM(); save_idxs=[X]).u[1]) == 1
    @test length(solve(jprob, SSAStepper(); save_idxs=[X]).u[1]) == 1
end

#Tests solution indexing.
let 
    for sol in deepcopy(sols)
        # Get u values.
        @test Int64(sol[X][1]) == 1
        @test Int64(sol[Y][1]) == 0
        @test sol[X] == sol[model.X] == sol[:X]
        @test sol[X2] == sol[model.X2] == sol[:X2]
        @test sol[[X,Y]] == sol[[model.X,model.Y]] == sol[[:X,:Y]]
        @test sol[[X2,Y]] == sol[[model.X2,model.Y]] == sol[[:X2,:Y]]
        @test sol[(X,Y)] == sol[(model.X,model.Y)] == sol[(:X,:Y)]   
        @test sol[(X2,Y)] == sol[(model.X2,model.Y)] == sol[(:X2,:Y)]   
        @test getu(sol, X)(sol) == getu(sol, model.X)(sol) == getu(sol, :X)(sol)    
        @test getu(sol, X2)(sol) == getu(sol, model.X2)(sol) == getu(sol, :X2)(sol)          

        # Get p values.
        @test getp(sol, kp)(sol) == getp(sol, model.kp)(sol) == getp(sol, :kp)(sol) == 1.0
        @test sol.ps[kp] == sol.ps[model.kp] == sol.ps[:kp] == 1.0    
    end
end

# Tests plotting.
let 
    for sol in deepcopy(sols[1:3])
        # Single variable.
        @test length(plot(sol; idxs = X).series_list) == 1
        @test length(plot(sol; idxs = X2).series_list) == 1
        @test length(plot(sol; idxs = model.X).series_list) == 1
        @test length(plot(sol; idxs = model.X2).series_list) == 1
        @test length(plot(sol; idxs = :X).series_list) == 1
        @test length(plot(sol; idxs = :X2).series_list) == 1

        # As vector.
        @test length(plot(sol; idxs = [X,Y]).series_list) == 2
        @test length(plot(sol; idxs = [X2,Y]).series_list) == 2
        @test length(plot(sol; idxs = [model.X,model.Y]).series_list) == 2
        @test length(plot(sol; idxs = [model.X2,model.Y]).series_list) == 2
        @test length(plot(sol; idxs = [:X,:Y]).series_list) == 2
        @test length(plot(sol; idxs = [:X2,:Y]).series_list) == 2

        # As tuple.
        @test length(plot(sol; idxs = (X, Y)).series_list) == 1
        @test length(plot(sol; idxs = (X2, Y)).series_list) == 1
        @test length(plot(sol; idxs = (model.X, model.Y)).series_list) == 1
        @test length(plot(sol; idxs = (model.X2, model.Y)).series_list) == 1
        @test length(plot(sol; idxs = (:X, :Y)).series_list) == 1
        @test length(plot(sol; idxs = (:X2, :Y)).series_list) == 1
    end    
end

# Tests solving for various inputs types.
let 
    u0_vals = [X => 1, Y => 0]
    tspan = (0.0, 10.0)
    p_vals = [kp => 1.0, kd => 0.2, k1 => 1.0, k2 => 2.0]

    u0_vals_2 = [model.X => 1, model.Y => 0]
    u0_vals_3 = [:X => 1, :Y => 0]
    p_vals_2 = [model.kp => 1.0, model.kd => 0.2, model.k1 => 1.0, model.k2 => 2.0]
    p_vals_3 = [:kp => 1.0, :kd => 0.2, :k1 => 1.0, :k2 => 2.0]

    oprob_2 = ODEProblem(model, u0_vals_2, tspan, p_vals_2)
    oprob_3 = ODEProblem(model, u0_vals_3, tspan, p_vals_3)
    sprob_2 = SDEProblem(model,u0_vals_2, tspan, p_vals_2)
    sprob_3 = SDEProblem(model,u0_vals_3, tspan, p_vals_3)
    dprob_2 = DiscreteProblem(model, u0_vals_2, tspan, p_vals_2)
    dprob_3 = DiscreteProblem(model, u0_vals_3, tspan, p_vals_3)
    jprob_2 = JumpProblem(model, dprob_2, Direct())
    jprob_3 = JumpProblem(model, dprob_3, Direct())
    nprob_2 = NonlinearProblem(model, u0_vals_2, p_vals_2)
    nprob_3 = NonlinearProblem(model, u0_vals_3, p_vals_3)
    
    @test solve(oprob, Tsit5()) == solve(oprob_2, Tsit5()) == solve(oprob_3, Tsit5())
    @test solve(sprob, ImplicitEM(); seed=1234) == solve(sprob_2, ImplicitEM(); seed=1234) == solve(sprob_3, ImplicitEM(); seed=1234)
    @test solve(jprob, SSAStepper(); seed=1234) == solve(jprob_2, SSAStepper(); seed=1234) == solve(jprob_3, SSAStepper(); seed=1234)
    @test solve(nprob, NewtonRaphson()) == solve(nprob_2, NewtonRaphson()) == solve(nprob_3, NewtonRaphson())
end