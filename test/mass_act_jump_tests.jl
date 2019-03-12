dotestmean   = true
doprintmeans = false
reltol       = .01          # required test accuracy
algs      = (Direct(),SortingDirect())

# run the given number of SSAs and return the mean
function runSSAs(jump_prob, Nsims, idx)
    Psamp = zeros(Int, Nsims)
    for i in 1:Nsims
        sol = solve(jump_prob, SSAStepper())
        Psamp[i] = sol[idx,end]
    end
    mean(Psamp)
end

function execute_test(u0, tf, rates, rs, Nsims, expected_avg, idx, test_name)
    prob = DiscreteProblem(u0, (0.0, tf), rates)

    for method in algs
        jump_prob = JumpProblem(prob, method, rs)
        avg_val = runSSAs(jump_prob, Nsims, idx)

        if dotestmean
            if doprintmeans
                println(test_name, ", method = ", typeof(method), ", mean = ", avg_val, ", act_mean = ", expected_avg)
            end
            @test abs(avg_val - expected_avg) < reltol * expected_avg
        end
    end

end

# nonlinear reaction test
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
execute_test(u0, tf, rates, rs, Nsims, expected_avg, 1, "Nonlinear rx test")


# DNA repression model
rs = @reaction_network ptype begin
    k1, DNA --> mRNA + DNA
    k2, mRNA --> mRNA + P
    k3, mRNA --> 0
    k4, P --> 0
    k5, DNA + P --> DNAR
    k6, DNAR --> DNA + P
end k1 k2 k3 k4 k5 k6
Nsims        = 8000
tf           = 1000.0
u0           = [1,0,0,0]
expected_avg = 5.926553750000000e+02
rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
execute_test(u0, tf, rates, rs, Nsims, expected_avg, 3, "DNA test")


# DNA repression model, mix of jump types
rs = @reaction_network ptype begin
    k1*DNA, 0 --> mRNA
    k2*mRNA, 0 --> P
    k3, mRNA --> 0
    k4, P --> 0
    k5, DNA + P --> DNAR
    k6, DNAR --> DNA + P
end k1 k2 k3 k4 k5 k6
Nsims        = 8000
tf           = 1000.0
u0           = [0,0,0,0]
u0[something(findfirst(isequal(:DNA),rs.syms),0)] = 1
expected_avg = 5.926553750000000e+02
rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
execute_test(u0, tf, rates, rs, Nsims, expected_avg, something(findfirst(isequal(:P),rs.syms), 0), "DNA mixed jump test")

# simple constant production with degratation
rs = @reaction_network pdtype begin
    1000., 0 --> A
    10, A --> 0
end
rates = nothing
Nsims        = 16000
tf           = 1.0
u0           = [0]
expected_avg = 1000.0 /10*(1. - exp(-10*tf))
execute_test(u0, tf, rates, rs, Nsims, expected_avg, 1, "Zero order test")


# this is just a test to make sure a mixture of jumps actually runs without crashes
network = @reaction_network rnType  begin
    0.01, (X,Y,Z) --> 0
    hill(X,3.,100.,-4), 0 --> Y
    hill(Y,3.,100.,-4), 0 --> Z
    hill(Z,4.5,100.,-4), 0 --> X
    hill(X,2.,100.,6), 0 --> R
    hill(Y,15.,100.,4)*0.002, R --> 0
    20, 0 --> S
    R*0.005, S --> SP
    0.01, SP + SP --> SP2
    0.05, SP2 --> 0
end;
prob = DiscreteProblem([200.,60.,120.,100.,50.,50.,50.], (0.,4000.))
for method in algs
    jump_prob = JumpProblem(prob, method, network)
    sol = solve(jump_prob,SSAStepper());
end

# make sure problem instantiation works when mass action jumps
# rate consts are a nontrivial expression on the params
network = @reaction_network begin
    1, X --> 2*X
    1/K, 2X --> X
end K
p = [1000.]
prob = DiscreteProblem([500.], (0.,100.), p)
jump_prob = JumpProblem(prob, Direct(), network)


# same as above
network = @reaction_network begin
    p1*p2, X --> Y
end p1 p2
p = [1.,2.]
prob = DiscreteProblem([10.,10.],(0.,10.), p)
jump_prob = JumpProblem(prob, Direct(), network)
