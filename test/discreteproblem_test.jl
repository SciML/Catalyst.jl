dotestmean   = true
doprintmeans = false
reltol       = .01          # required test accuracy
algs      = (Direct(),)

# run the given number of SSAs and return the mean
function runSSAs(jump_prob, Nsims, idx)
    Psamp = zeros(Int, Nsims)
    for i in 1:Nsims
        sol = solve(jump_prob, SSAStepper())
        Psamp[i] = sol[idx,end]
    end
    mean(Psamp)
end

function execute_test(prob, u0, tf, rates, rs, Nsims, expected_avg, idx, test_name)
    for method in algs
        jump_prob = JumpProblem(prob, method, rs, save_positions=(false,false))
        avg_val = runSSAs(jump_prob, Nsims, idx)

        if dotestmean
            if doprintmeans
                println(test_name, ", method = ", typeof(method), ", mean = ", avg_val, ", act_mean = ", expected_avg)
            end
            @test abs(avg_val - expected_avg) < reltol * expected_avg
        end
    end

end

# full macro
Nsims = 32000
tf = .01
u0 = [200, 100, 150]
expected_avg = 84.876015624999994
rs = @reaction_network dtype begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
rates = [1., 2., .5, .75, .25]
prob = DiscreteProblem(u0, (0.0, tf), rates)
execute_test(prob, u0, tf, rates, rs, Nsims, expected_avg, 1, "Nonlinear rx test (no syms)")

prob = DiscreteProblem(rs, u0, (0.0, tf), rates)
execute_test(prob, u0, tf, rates, rs, Nsims, expected_avg, 1, "Nonlinear rx test (syms)")

# min_network macro
rs = @min_reaction_network dtype begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
addjumps!(rs)
prob = DiscreteProblem(u0, (0.0, tf), rates)
execute_test(prob, u0, tf, rates, rs, Nsims, expected_avg, 1, "Nonlinear rx test, min network (no syms)")

prob = DiscreteProblem(rs, u0, (0.0, tf), rates)
execute_test(prob, u0, tf, rates, rs, Nsims, expected_avg, 1, "Nonlinear rx test, min network (syms)")
