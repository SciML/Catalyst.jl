### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, LinearAlgebra, NonlinearSolve, OrdinaryDiffEqTsit5, OrdinaryDiffEqDefault, OrdinaryDiffEqVerner, SteadyStateDiffEq, StochasticDiffEq, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(123456)
seed = rand(rng, 1:100)

# Fetch test networks.
include("../test_networks.jl")

### Basic Tests ###

# Tests basic functionality on system with known conservation laws.
let
    rn = @reaction_network begin
        (1, 2), A + B <--> C
        (3, 2), D <--> E
        (0.1, 0.2), E <--> F
        6, F --> G
        7, H --> G
        5, 0 --> K
        3, A + B --> Z + C
    end

    S = netstoichmat(rn)
    C = conservationlaws(S)
    @test size(C, 1) == 3
    b = [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
    @test any(b == C[i, :] for i in 1:size(C, 1))

    # For the A + B <--> C subsystem one of these must occur
    # as a conservation law.
    D = [1 -1 0 0 0 0 0 0 0 0;
         -1 1 0 0 0 0 0 0 0 0
         0 1 1 0 0 0 0 0 0 0]
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))

    C = conservationlaws(rn)
    @test size(C, 1) == 3
    @test Catalyst.get_networkproperties(rn).nullity == 3
    @test any(b == C[i, :] for i in 1:size(C, 1))
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))
end

# Tests conservation law computation on large number of networks where we know which have conservation laws.
let
    # networks for which we know there is no conservation laws.
    Cs_standard = map(conservationlaws, reaction_networks_standard)
    Cs_hill = map(conservationlaws, reaction_networks_hill)
    @test all(size(C, 1) == 0 for C in Cs_standard)
    @test all(size(C, 1) == 0 for C in Cs_hill)

    # Networks for which there are known conservation laws (stored in `reaction_network_conslaws`).
    function consequiv(A, B)
        rank([A; B]) == rank(A) == rank(B)
    end
    Cs_constraint = map(conservationlaws, reaction_networks_conserved)
    @test all(consequiv.(Matrix{Int}.(Cs_constraint), reaction_network_conslaws))
end

# Tests additional conservation law-related functions.
let
    rn = @reaction_network begin
        (k1, k2), X1 <--> X2
        (k3, k4), X3 <--> X4
    end
    cons_laws = conservationlaws(rn)
    cons_eqs = conservedequations(rn)
    cons_laws_constants = conservationlaw_constants(rn)
    conserved_quantity = conservedquantities(cons_laws[1, :], unknowns(rn)[1])

    @test sum(cons_laws) == 4
    @test size(cons_laws) == (2, 4)
    @test length(cons_eqs) == 2
    @test length(conserved_quantity) == 4
    @test length(cons_laws_constants) == 2
    @test count(isequal.(conserved_quantity, Num(0))) == 2
end

# Tests that `conservationlaws`'s caches something.
let
    # Creates network with/without cached conservation laws.
    rn = @reaction_network rn begin
        (k1,k2), X1 <--> X2
    end
    rn_cached = deepcopy(rn)
    conservationlaws(rn_cached)

    # Checks that equality is correct (currently equality does not consider network property caching).
    @test rn_cached == rn
    @test Catalyst.get_networkproperties(rn_cached) != Catalyst.get_networkproperties(rn)
end

### Simulation & Solving Tests ###

# Test conservation law elimination.
let
    # Declares the model
    rn = @reaction_network begin
        (k1, k2), A + B <--> C
        (m1, m2), D <--> E
        b12, F1 --> F2
        b23, F2 --> F3
        b31, F3 --> F1
    end
    @unpack A, B, C, D, E, F1, F2, F3, k1, k2, m1, m2, b12, b23, b31 = rn
    sps = species(rn)
    u0 = [A => 10.0, B => 10.0, C => 0.0, D => 10.0, E => 0.0, F1 => 8.0, F2 => 0.0,
        F3 => 0.0]
    p = [k1 => 1.0, k2 => 0.1, m1 => 1.0, m2 => 2.0, b12 => 1.0, b23 => 2.0, b31 => 0.1]
    tspan = (0.0, 20.0)

    # Simulates model using ODEs and checks that simulations are identical.
    osys = complete(convert(ODESystem, rn; remove_conserved = true))
    oprob1 = ODEProblem(osys, u0, tspan, p)
    oprob2 = ODEProblem(rn, u0, tspan, p)
    oprob3 = ODEProblem(rn, u0, tspan, p; remove_conserved = true)
    osol1 = solve(oprob1, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    osol2 = solve(oprob2, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    osol3 = solve(oprob3, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    @test osol1[sps] ≈ osol2[sps] ≈ osol3[sps]

    # Checks that steady states found using nonlinear solving and steady state simulations are identical.
    nsys = convert(NonlinearSystem, rn; remove_conserved = true, conseqs_remake_warn = false)
    nsys_ss = structural_simplify(nsys)
    nsys = complete(nsys)
    nprob1 = NonlinearProblem{true}(nsys, u0, p)
    nprob1b = NonlinearProblem{true}(nsys_ss, u0, p)
    nprob2 = NonlinearProblem(rn, u0, p; remove_conserved = true, conseqs_remake_warn = false)
    nprob2b = NonlinearProblem(rn, u0, p; remove_conserved = true, conseqs_remake_warn = false, 
        structural_simplify = true)
    nsol1 = solve(nprob1, NewtonRaphson(); abstol = 1e-8)
    nsol1b = solve(nprob1b, NewtonRaphson(); abstol = 1e-8)
    nsol2 = solve(nprob2, NewtonRaphson(); abstol = 1e-8)
    nsol2b = solve(nprob2b, NewtonRaphson(); abstol = 1e-8)
    # Nonlinear problems cannot find steady states properly without removing conserved species.

    ssprob1 = SteadyStateProblem{true}(osys, u0, p)
    ssprob2 = SteadyStateProblem(rn, u0, p)
    ssprob3 = SteadyStateProblem(rn, u0, p; remove_conserved = true)
    sssol1 = solve(ssprob1, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    sssol2 = solve(ssprob2, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    sssol3 = solve(ssprob3, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    @test nsol1[sps] ≈ nsol1b[sps] 
    @test nsol1[sps] ≈ nsol2[sps] 
    @test nsol1[sps] ≈ nsol2b[sps] 
    @test nsol1[sps] ≈ sssol1[sps]
    @test nsol1[sps] ≈ sssol2[sps] 
    @test nsol1[sps] ≈ sssol3[sps]

    # Creates SDEProblems using various approaches.
    u0_sde = [A => 100.0, B => 20.0, C => 5.0, D => 10.0, E => 3.0, F1 => 8.0, F2 => 2.0,
        F3 => 20.0]
    ssys = complete(convert(SDESystem, rn; remove_conserved = true))
    sprob1 = SDEProblem(ssys, u0_sde, tspan, p)
    sprob2 = SDEProblem(rn, u0_sde, tspan, p)
    sprob3 = SDEProblem(rn, u0_sde, tspan, p; remove_conserved = true)

    # Checks that the SDEs f and g function evaluates to the same thing.
    ind_us = ModelingToolkit.get_unknowns(ssys)
    us = ModelingToolkit.get_unknowns(rn)
    ind_uidxs = findall(in(ind_us), us)
    u1 = copy(sprob1.u0)
    u2 = sprob2.u0
    u3 = copy(sprob3.u0)
    du1 = similar(u1)
    du2 = similar(u2)
    du3 = similar(u3)
    g1 = zeros(length(u1), numreactions(rn))
    g2 = zeros(length(u2), numreactions(rn))
    g3 = zeros(length(u3), numreactions(rn))
    sprob1.f(du1, sprob1.u0, sprob1.p, 1.0)
    sprob2.f(du2, sprob2.u0, sprob2.p, 1.0)
    sprob3.f(du3, sprob3.u0, sprob3.p, 1.0)
    @test du1 ≈ du2[ind_uidxs] ≈ du3
    sprob1.g(g1, sprob1.u0, sprob1.p, 1.0)
    sprob2.g(g2, sprob2.u0, sprob2.p, 1.0)
    sprob3.g(g3, sprob3.u0, sprob3.p, 1.0)
    @test g1 ≈ g2[ind_uidxs, :] ≈ g3
end

# Tests simulations for various input types (using X, rn.X, and :X forms).
# Tests that conservation laws can be generated for system with non-default parameter types.
let
    # Prepares the model.
    rn = @reaction_network rn begin
        @parameters kB::Int64
        (kB,kD), X + Y <--> XY
    end
    sps = species(rn)
    @unpack kB, kD, X, Y, XY = rn

    # Creates `ODEProblem`s using three types of inputs. Checks that solutions are identical.
    u0_1 = [X => 2.0, Y => 3.0, XY => 4.0]
    u0_2 = [rn.X => 2.0, rn.Y => 3.0, rn.XY => 4.0]
    u0_3 = [:X => 2.0, :Y => 3.0, :XY => 4.0]
    ps = (kB => 2, kD => 1.5)
    oprob1 = ODEProblem(rn, u0_1, 10.0, ps; remove_conserved = true)
    oprob2 = ODEProblem(rn, u0_2, 10.0, ps; remove_conserved = true)
    oprob3 = ODEProblem(rn, u0_3, 10.0, ps; remove_conserved = true)
    @test solve(oprob1)[sps] ≈ solve(oprob2)[sps] ≈ solve(oprob3)[sps]
end

# Tests conservation laws in SDE simulation.
let
    # Creates `SDEProblem`s.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    u0 = Dict([:X1 => 100.0, :X2 => 120.0])
    ps = [:k1 => 0.2, :k2 => 0.15]
    sprob = SDEProblem(rn, u0, 10.0, ps; remove_conserved = true)

    # Checks that conservation laws hold in all simulations.
    sol = solve(sprob, ImplicitEM(); seed)
    @test sol[:X1] + sol[:X2] ≈ sol[rn.X1 + rn.X2] ≈ fill(u0[:X1] + u0[:X2], length(sol.t))
end

# Checks that the conservation law parameter's value can be changed in simulations.
let
    # Prepares `ODEProblem`s.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    osys = complete(convert(ODESystem, rn; remove_conserved = true))
    u0_1 = [osys.X1 => 1.0, osys.X2 => 1.0]
    u0_2 = [osys.X1 => 1.0]
    ps_1 = [osys.k1 => 2.0, osys.k2 => 3.0]
    ps_2 = [osys.k1 => 2.0, osys.k2 => 3.0, osys.Γ => [4.0]]
    oprob1 = ODEProblem(osys, u0_1, 10.0, ps_1)
    oprob2 = ODEProblem(osys, u0_2, 10.0, ps_2)

    # Checks that the solutions generates the correct conserved quantities.
    sol1 = solve(oprob1; saveat = 1.0)
    sol2 = solve(oprob2; saveat = 1.0)
    @test all(sol1[osys.X1 + osys.X2] .== 2.0)
    @test all(sol2[osys.X1 + osys.X2] .== 4.0)
end

### Problem `remake`ing and Updating Tests ###

# Tests system problem updating when conservation laws are eliminated.
# Checks that the correct values are used after the conservation law species are updated.
# Here is an issue related to the broken tests: https://github.com/SciML/Catalyst.jl/issues/952
let
    # Create model and fetch the conservation parameter (Γ).
    t = default_t()
    @parameters k1 k2
    @species X1(t) X2(t)
    rxs = [
        Reaction(k1, [X1], [X2]),
        Reaction(k2, [X2], [X1])
    ]
    @named rs = ReactionSystem(rxs, t)
    osys = complete(convert(ODESystem, complete(rs); remove_conserved = true))
    @unpack Γ = osys

    # Creates the various problem types.
    u0 = [X1 => 1.0, X2 => 2.0]
    ps = [k1 => 0.1, k2 => 0.2]
    oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)

    # Check `ODEProblem` content.
    @test oprob[X1] == 1.0
    @test oprob[X2] == 2.0
    @test oprob.ps[k1] == 0.1
    @test oprob.ps[k2] == 0.2
    @test oprob.ps[Γ[1]] == 3.0

    # Update problem parameters using `remake`.
    oprob_new = remake(oprob; p = [k1 => 0.3, k2 => 0.4])
    @test oprob_new.ps[k1] == 0.3
    @test oprob_new.ps[k2] == 0.4
    integrator = init(oprob_new, Tsit5())
    @test integrator.ps[k1] == 0.3
    @test integrator.ps[k2] == 0.4

    # Update problem parameters using direct indexing.
    oprob.ps[k1] = 0.5
    oprob.ps[k2] = 0.6
    @test oprob.ps[k1] == 0.5
    @test oprob.ps[k2] == 0.6
    integrator = init(oprob, Tsit5())
    @test integrator.ps[k1] == 0.5
    @test integrator.ps[k2] == 0.6
end

### Jacobian Tests ###

# Checks that conservation law elimination generates a system with a non-singular Jacobian.
# Does this for different system types (ODE, SDE, and Nonlinear).
# Checks singularity by checking whether the Jacobian have high enough a condition number 
# (when evaluated at the steady state).
let
    # Creates the model (contains both conserved and non-conserved species).
    rn = @reaction_network begin
        (p,d), 0 <--> X
        d, XY --> Y
        (k1,k2), X + Y <--> XY
        (k3,k4), 2Y <--> Y2
    end

    # Finds a steady state at which we will compute the Jacobian.
    u0_ss_init = [:X => 0.0, :Y => 4.0, :XY => 0.0, :Y2 => 0.0]
    ps = [:p => 1.0, :d => 0.5, :k1 => 1.0, :k2 => 2.0, :k3 => 3.0, :k4 => 4.0]
    ssprob = SteadyStateProblem(rn, u0_ss_init, ps)
    ss = solve(ssprob, DynamicSS(Vern7()))
    ss = [sp => ss[sp] for sp in unknowns(rn)]

    # Creates a function which evaluates whether the Jacobian is singular.
    # Singularity means infinite condition number (here it is about 1e17).
    function is_singular(prob; infthres = 1e12)
        J = zeros(length(prob.u0), length(prob.u0))
        ModelingToolkit.is_time_dependent(prob) ? prob.f.jac(J, prob.u0, prob.p, 0.0) : prob.f.jac(J, prob.u0, prob.p)
        return cond(J) > infthres
    end

    # Creates problems (with Jacobian and W/O conservation law elimination) which singularity we wish to test.
    oprob = ODEProblem(rn, ss, 1.0, ps; jac = true)
    oprob_rc = ODEProblem(rn, ss, 1.0, ps; jac = true, remove_conserved = true)
    sprob = SDEProblem(rn, ss, 1.0, ps; jac = true)
    sprob_rc = SDEProblem(rn, ss, 1.0, ps; jac = true, remove_conserved = true)
    nlprob = NonlinearProblem(rn, ss, ps; jac = true)
    nlprob_rc = NonlinearProblem(rn, ss, ps; jac = true, remove_conserved = true, conseqs_remake_warn = false)

    # Checks that removing conservation laws generates non-singular Jacobian (and else that it is singular).
    @test is_singular(oprob) == true
    @test !is_singular(oprob_rc) == true
    @test is_singular(sprob) == true
    @test !is_singular(sprob_rc) == true
    @test is_singular(nlprob) == true
    @test !is_singular(nlprob_rc) == true
end

# Goes through a chain of updating of conservation law constants/species, checking that
# the values of relevant quantities are correct after each step.
# Generally, if `Γ` has not been explicitly updated, it will be updated to acomodate new species
# values. If it has been explicitly updated, the corresponding eliminated quantity will have its
# value updated to acomodate new Γ/species values (however, we have to manually designate this by setting it to `nothing`).
# Also checks that quantities are correctly updated in integrators and solutions derived from problems.
let
    # Prepares the problem inputs and computes the conservation equation.
    rn = @reaction_network begin
        (k1,k2), 2X1 <--> X2
        (k3,k4), X1 + X2 <--> X3
    end
    @unpack X1, X2, X3 = rn
    u0 = [X1 => 1.0, X2 => 1.0, X3 => 1.0]
    ps = [:k1 => 0.1, :k2 => 0.2, :k3 => 0.3, :k4 => 0.4]
    conserved_quantity = conservationlaw_constants(rn)[1].rhs

    # Loops through the tests for different problem types.
    oprob = ODEProblem(rn, u0, 1.0, ps; remove_conserved = true)
    sprob = SDEProblem(rn, u0, 1.0, ps; remove_conserved = true)
    # nlprob = NonlinearProblem(rn, u0, ps; remove_conserved = true, conseqs_remake_warn = false)
    @test_broken false # Cases of `remake` does not work for `NonlinearSystem`s with `remove_conserved = true`.
    for (prob_old, solver) in zip([oprob, sprob], [Tsit5(), ImplicitEM()])
        # For a couple of iterations, updates the problem, ensuring that when a species is updated:
        # - Only that species and the conservation constant have their values updated.
        # The `≈` is because sometimes the computed values will not be fully exact.
        for _ = 1:3
            # Updates X2, checks the values of all species and Γ, then resets `prob_old`.
            X2_new = rand(rng, 1.0:10.0)
            prob_new = remake(prob_old; u0 = [:X2 => X2_new])
            @test prob_old[:X1] ≈ prob_new[:X1]
            @test X2_new ≈ prob_new[:X2]
            @test prob_old[:X3] ≈ prob_new[:X3]
            @test substitute(conserved_quantity, Dict([X1 => prob_old[X1], X2 => X2_new, X3 => prob_old[X3]])) ≈ prob_new.ps[:Γ][1]
            prob_old = prob_new

            # Updates X3, checks the values of all species and Γ, then resets `prob_old`.
            X3_new = rand(rng, 1.0:10.0)
            prob_new = remake(prob_old; u0 = [:X3 => X3_new])
            @test prob_old[:X1] ≈ prob_new[:X1]
            @test prob_old[:X2] ≈ prob_new[:X2]
            @test X3_new ≈ prob_new[:X3]
            @test substitute(conserved_quantity, Dict([X1 => prob_old[X1], X2 => prob_old[X2], X3 => X3_new])) ≈ prob_new.ps[:Γ][1]
            prob_old = prob_new

            # Checks that integrator and solutions have identical content to problem.
            integrator = init(prob_new, solver)
            sol = solve(prob_new, solver; maxiters = 1, verbose = false)
            @test prob_new.ps[:Γ][1] == integrator.ps[:Γ][1] == sol.ps[:Γ][1]
            if ModelingToolkit.is_time_dependent(prob_new)
                @test prob_new[:X1] == integrator[:X1] == sol[:X1][1]
                @test prob_new[:X2] == integrator[:X2] == sol[:X2][1]
                @test prob_new[:X3] == integrator[:X3] == sol[:X3][1]
            else
                @test prob_new[:X1] == integrator[:X1]
                @test prob_new[:X2] == integrator[:X2]
                @test prob_new[:X3] == integrator[:X3]
            end
        end

        # Similarly, but now also updates the conservation constant. Here, once Γ has been updated:
        # - The conservation law constant will be kept fixed, and secondary updates are made to the
        # eliminated species.
        # Assumes that X3 is the eliminated species. In most updates, designate its as `nothing` (to
        # ensure that its value is updated to acommodate the new conservation law).
        # The random Γ is ensured to be large enough not to generate negative values in the eliminated species.
        for _ in 1:3
            # Updates Γ, checks the values of all species and Γ, then resets `prob_old`.
            Γ_new = substitute(conserved_quantity, Dict([X1 => prob_old[X1], X2 => prob_old[X2], X3 => 0])) + rand(rng, 0.0:5.0)
            prob_new = remake(prob_old; u0 = [:X3 => nothing], p = [:Γ => [Γ_new]])
            @test prob_old[:X1] ≈ prob_new[:X1]
            @test prob_old[:X2] ≈ prob_new[:X2]
            @test Γ_new ≈ prob_new.ps[:Γ][1]
            @test substitute(conserved_quantity, Dict([X1 => prob_old[X1], X2 => prob_old[X2], X3 => prob_new[X3]])) ≈ prob_new.ps[:Γ][1]
            prob_old = prob_new

            # Updates X1 (non-eliminated species), checks the values of all species and Γ, then resets `prob_old`.
            # Note that now, `X3` will have its value modified (and `Γ` remains unchanged).
            X1_new = rand(rng, 1.0:10.0)
            prob_new = remake(prob_old; u0 = [:X1 => X1_new, :X3 => nothing])
            @test X1_new ≈ prob_new[:X1]
            @test prob_old[:X2] ≈ prob_new[:X2]
            @test prob_old.ps[:Γ][1] ≈ prob_new.ps[:Γ][1]
            @test substitute(conserved_quantity, Dict([X1 => X1_new, X2 => prob_old[X2], X3 => prob_new[X3]])) ≈ prob_new.ps[:Γ][1]
            prob_old = prob_new

            # Updates X3 (the eliminated species). Since we reset Γ, this has its value modified to 
            # accommodate the new value of X3.
            X3_new = rand(rng, 1.0:10.0)
            prob_new = remake(prob_old; u0 = [:X3 => X3_new], p = [:Γ => nothing])
            @test prob_old[:X1] ≈ prob_new[:X1]
            @test prob_old[:X2] ≈ prob_new[:X2]
            @test X3_new ≈ prob_new[:X3]
            @test substitute(conserved_quantity, Dict([X1 => prob_old[X1], X2 => prob_old[X2], X3 => X3_new])) ≈ prob_new.ps[:Γ][1]
            prob_old = prob_new

            # Checks that integrator and solutions have identical content to problem.
            integrator = init(prob_new, solver)
            sol = solve(prob_new, solver; maxiters = 1, verbose = false)
            @test prob_new.ps[:Γ][1] == integrator.ps[:Γ][1] == sol.ps[:Γ][1]
            if ModelingToolkit.is_time_dependent(prob_new)
                @test prob_new[:X1] == integrator[:X1] == sol[:X1][1]
                @test prob_new[:X2] == integrator[:X2] == sol[:X2][1]
                @test prob_new[:X3] == integrator[:X3] == sol[:X3][1]
            else
                @test prob_new[:X1] == integrator[:X1]
                @test prob_new[:X2] == integrator[:X2]
                @test prob_new[:X3] == integrator[:X3]
            end
        end
    end
end

### Other Tests ###

# Checks that `JumpSystem`s with conservation laws cannot be generated.
let
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    @test_throws ArgumentError convert(JumpSystem, rn; remove_conserved = true)
end

# Checks that `conserved` metadata is added correctly to parameters.
# Checks that the `isconserved` getter function works correctly.
let
    # Creates ODESystem with conserved quantities.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k1,k2), Y1 <--> Y2
    end
    osys = convert(ODESystem, rs; remove_conserved = true)

    # Checks that the correct parameters have the `conserved` metadata.
    @test Catalyst.isconserved(osys.Γ[1])
    @test Catalyst.isconserved(osys.Γ[2])
    @test !Catalyst.isconserved(osys.k1)
    @test !Catalyst.isconserved(osys.k2)
end

# Conservation law simulations for vectorised species.
let
    # Prepares the model.
    t = default_t()
    @species (X(t))[1:2]
    @parameters k[1:2]
    rxs = [
        Reaction(k[1], [X[1]], [X[2]]),
        Reaction(k[2], [X[2]], [X[1]])
    ]
    @named rs = ReactionSystem(rxs, t)
    rs = complete(rs)

    # Checks that simulation reaches known equilibrium.
    @test_broken false # Currently broken on MTK .
    # u0 = [:X => [3.0, 9.0]]
    # ps = [:k => [1.0, 2.0]]
    # oprob = ODEProblem(rs, u0, (0.0, 1000.0), ps; remove_conserved = true)
    # sol = solve(oprob, Vern7())
    # @test sol[X[1]][end] ≈ 8.0
    # @test sol[X[2]][end] ≈ 4.0
end

# Check conservation law elimination warnings (and the disabling of these) for `NonlinearSystem`s
# and `NonlinearProblem`s.
let
    # Create models.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k3,k4), X1 + X2 <--> 2X3
    end
    u0 = [:X1 => 1.0, :X2 => 2.0, :X3 => 3.0]
    ps = [:k1 => 0.1, :k2 => 0.2, :k3 => 0.3, :k4 => 0.4]
    
    # Checks that the warning si given and can be supressed for the variosu cases.
    @test_nowarn convert(NonlinearSystem, rn; remove_conserved = true, conseqs_remake_warn = false)
    @test_logs (:warn, r"Note, when constructing*") convert(NonlinearSystem, rn; remove_conserved = true, conseqs_remake_warn = true)
    @test_nowarn NonlinearProblem(rn, u0, ps; remove_conserved = true, conseqs_remake_warn = false)
    @test_logs (:warn, Catalyst.NONLIN_PROB_REMAKE_WARNING) NonlinearProblem(rn, u0, ps; remove_conserved = true, conseqs_remake_warn = true)

    # @variables X1 X2 X3
    # @parameters k1 k2 k3 k4 Γ[1:1] = missing [guess = ones(1)]
    # eqs = [
    #     0 ~ -k1*X1 + k2*X2 - k3*X1*X2 + (1//2)*k4*((-X1 - X2 + Γ[1])^2),
    #     0 ~ k1*X1 - k2*X2 - k3*X1*X2 + (1//2)*k4*((-X1 - X2 + Γ[1])^2),
    #     0 ~ -X1 - X2 - X3 + Γ[1]
    # ]
    # initeqs = [Γ[1] ~ Initial(X1) + Initial(X3) + Initial(X2)]
    # @named nlsys = NonlinearSystem(eqs, [X1, X2, X3], [k1, k2, k3, k4, Γ]; 
    #     initialization_eqs = initeqs)

    
    # WITHOUT structural_simplify
    nlsys = convert(NonlinearSystem, rn; remove_conserved = true, 
        conseqs_remake_warn = false)
    nlsys1 = complete(nlsys)
    nlprob1 = NonlinearProblem(nlsys1, u0, ps)

    @test nlprob1[:X1] == 1.0
    @test nlprob1[:X2] == 2.0
    @test nlprob1[:X3] == 3.0
    @test nlprob1.ps[:Γ][1] == 6.0
    integ1 = init(nlprob1, NewtonRaphson())
    @test integ1[:X1] == 1.0
    @test integ1[:X2] == 2.0
    @test integ1[:X3] == 3.0
    @test integ1.ps[:Γ][1] == 6.0

    nlprob1b = remake(nlprob1; u0 = [:X3 => nothing], p = [:Γ => [4.0]])
    @test nlprob1b[:X1] == 1.0   
    @test nlprob1b[:X2] == 2.0   
    @test_broken nlprob1b[:X3] == 1.0
    @test nlprob1b.ps[:Γ][1] == 4.0
    integ1 = init(nlprob1b, NewtonRaphson())
    @test integ1[:X1] == 1.0
    @test integ1[:X2] == 2.0
    @test_broken integ1[:X3] == 1.0
    @test integ1.ps[:Γ][1] == 4.0

    nlprob1c = remake(nlprob1; u0 = [:X2 => nothing], p = [:Γ => [4.0]])
    @test nlprob1c[:X1] == 1.0
    @test_broken nlprob1c[:X2] == 0.0
    @test nlprob1c[:X3] == 3.0
    @test nlprob1c.ps[:Γ][1] == 4.0
    integ1 = init(nlprob1b, NewtonRaphson())
    @test integ1[:X1] == 1.0
    @test_broken integ1[:X2] == 0.0
    @test integ1[:X3] == 3.0
    @test integ1.ps[:Γ][1] == 4.0
    
    # WITH structural_simplify
    nlsys2 = structural_simplify(nlsys)
    nlprob2 = NonlinearProblem(nlsys2, u0, ps)

    @test nlprob2[:X1] == 1.0
    @test nlprob2[:X2] == 2.0
    @test nlprob2[:X3] == 3.0
    @test nlprob2.ps[:Γ][1] == 6.0
    integ2 = init(nlprob2, NewtonRaphson())
    @test integ2[:X1] == 1.0
    @test integ2[:X2] == 2.0
    @test integ2[:X3] == 3.0
    @test integ2.ps[:Γ][1] == 6.0

    nlprob2b = remake(nlprob2; u0 = [:X3 => nothing], p = [:Γ => [4.0]])
    @test nlprob2b[:X1] == 1.0   
    @test nlprob2b[:X2] == 2.0   
    @test nlprob2b[:X3] == 1.0   
    @test nlprob2b.ps[:Γ][1] == 4.0 
    integ2 = init(nlprob2b, NewtonRaphson())
    @test integ2[:X1] == 1.0  
    @test integ2[:X2] == 2.0  
    @test integ2[:X3] == 1.0  
    @test integ2.ps[:Γ][1] == 4.0 

    nlprob2c = remake(nlprob2; u0 = [:X2 => nothing], p = [:Γ => [4.0]])
    @test nlprob2c[:X1] == 1.0   
    @test_broken nlprob2c[:X2] == 0.0 
    @test_broken nlprob2c[:X3] == 3.0
    @test nlprob2c.ps[:Γ][1] == 4.0 
    integ2 = init(nlprob2c, NewtonRaphson())
    @test integ2[:X1] == 1.0  
    @test_broken integ2[:X2] == 0.0 
    @test_broken integ2[:X3] == 3.0
    @test integ2.ps[:Γ][1] == 4.0 
end