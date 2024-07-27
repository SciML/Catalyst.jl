### ODE & SDE Assembly ###

"""
    oderatelaw(rx; combinatoric_ratelaw = true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated ODEs for the reaction. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X(t)^2*Y(t)*Z(t)`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * (X(t)^2/2) * (Y(t)^3/6)`.

Notes:
- Allocates
- `combinatoric_ratelaw = true` uses factorial scaling factors in calculating the
 rate law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
 `combinatoric_ratelaw = false` then the ratelaw is `k*S^2`, i.e. the scaling
 factor is ignored.
"""
function oderatelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate

    # if the stoichiometric coefficients are not integers error if asking to scale rates
    !all(s -> s isa Union{Integer, Symbolic}, substoich) &&
        (combinatoric_ratelaw == true) &&
        error("Non-integer stoichiometric coefficients require the combinatoric_ratelaw = false keyword to oderatelaw, or passing combinatoric_ratelaws = false to convert or ODEProblem.")

    # if we use mass action, update the rate law according to the reaction substrates
    if !only_use_rate
        # `coef` stores the numeric part (i.e. division by factorial of stoichiometries)
        # the symbolic part (i.e. multiplication by substrate concentrations) is appended directly 
        # to the rate law (`rl`)
        # cases with stoichiometries 1 are handled more efficiently using separate routines
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            combinatoric_ratelaw && (coef *= factorial(stoich))
            rl *= isequal(stoich, one(stoich)) ? substrates[i] : substrates[i]^stoich
        end
        combinatoric_ratelaw && (!isequal(coef, one(coef))) && (rl /= coef)
    end
    rl
end

# function returning `true` for species which shouldn't change from the reactions, 
# including non-species variables
drop_dynamics(sp) = isconstant(sp) || isbc(sp) || (!isspecies(sp))

# for a `ReactionSystem`, create a vector with all the right-hand side expressions in the
# corresponding ODE system (used by the `assemble_drift` function)
function assemble_oderhs(rs, isps; combinatoric_ratelaws = true, remove_conserved = false)
    # initiates the right-hand side equations as expressions with `0` only.
    species_to_idx = Dict(x => i for (i, x) in enumerate(isps))
    rhsvec = Any[0 for _ in isps]

    # if we have eliminated conservation laws, create a dictionary taking each species eliminated
    # in this process to the expression it is replaced by (if not, creates an empty dict)
    nps = get_networkproperties(rs)
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    # loops through all reactions, for each, updates all relevant right-hand side equations
    for rx in get_rxs(rs)
        # computes the reaction's rate late
        # if conserved quantities have been removed, replaces these in the rate law
        rl = oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
        remove_conserved && (rl = substitute(rl, depspec_submap))

        # for each species which participates in the reaction, update its equation
        for (spec, stoich) in rx.netstoich
            # dependent species don't get an ODE, so are skipped
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            # if the relevant RHS equation is 0 (probably because it has not been
            # modified yet), we use a different routine (compared to if it is non-zero)
            if _iszero(rhsvec[i])
                if stoich isa Symbolic
                    rhsvec[i] = stoich * rl
                else
                    # If the stoichiometry is +/- 1, initiate using the (potentially negative)
                    # rate law directly. If not, initiate with rate law * stoichiometry.
                    signedrl = (stoich > zero(stoich)) ? rl : -rl
                    rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
                end
            else
                if stoich isa Symbolic
                    rhsvec[i] += stoich * rl
                else
                    # If the stoichiometry is negative, we subtract the rate law (multiplied by the
                    # absolute value of the stoichiometry) from the RHS. If not, add it.
                    Δspec = isone(abs(stoich)) ? rl : abs(stoich) * rl
                    rhsvec[i] = (stoich > zero(stoich)) ? (rhsvec[i] + Δspec) :
                                (rhsvec[i] - Δspec)
                end
            end
        end
    end

    rhsvec
end

# For a `ReactionSystem`, creates a vector with all the differential equations corresponding to
# the ODEs generated by the reaction rate equation. The `as_odes` argument determine whether the LHS should
# be differential (when generating e.g. a `ODESystem`) or `0` (when generating e.g. a NonlinearSystem).
function assemble_drift(rs, ind_sps; combinatoric_ratelaws = true, as_odes = true,
                        include_zero_odes = true, remove_conserved = false)
    rhsvec = assemble_oderhs(rs, ind_sps; combinatoric_ratelaws, remove_conserved)
    if as_odes
        D = Differential(get_iv(rs))
        eqs = [Equation(D(x), rhs)
               for (x, rhs) in zip(ind_sps, rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
    else
        eqs = [Equation(0, rhs) for rhs in rhsvec if (include_zero_odes || (!_iszero(rhs)))]
    end
    eqs
end

# this doesn't work with coupled equations currently
# assemblies the diffusion equations terms used in the chemical Langevin equation SDE
function assemble_diffusion(rs, us, ind_sps; combinatoric_ratelaws = true,
                            remove_conserved = false)
    # as BC species should ultimately get an equation, we include them in the noise matrix
    num_bcus = count(isbc, get_unknowns(rs))

    # we make a matrix sized by the number of reactions
    eqs = Matrix{Any}(undef, length(sts) + num_bcsts, length(get_rxs(rs)))
    eqs .= 0
    species_to_idx = Dict((x => i for (i, x) in enumerate(ispcs)))

    # if we have eliminated conservation laws, create a dictionary taking each species eliminated
    # in this process to the expression it is replaced by (if not, creates an empty dict)
    nps = get_networkproperties(rs)
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    # loops through all reactions, for each, updates relevant noise terms in the noise matrix
    for (j, rx) in enumerate(get_rxs(rs))
        # Computes the reaction's (noise) rate law. This is the sqrt of the rate law, multiplied
        # by any noise scaling term. If conserved quantities have been removed, these are substituted in.
        rlsqrt = sqrt(abs(oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)))
        hasnoisescaling(rx) && (rlsqrt *= getnoisescaling(rx))
        remove_conserved && (rlsqrt = substitute(rlsqrt, depspec_submap))

        # for each species which participate in the reaction, updates its corresponding noise terms
        for (spec, stoich) in rx.netstoich
            # dependent species don't get an equation
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            # Different routine if the stoichiometry is a symbolic (e.g. a parameter)
            # if not it might e +/- 1, in which case we initiate the term differently
            if stoich isa Symbolic
                eqs[i, j] = stoich * rlsqrt
            else
                signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
                eqs[i, j] = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
            end
        end
    end
    eqs
end

### Jumps Assembly ###

"""
 jumpratelaw(rx; combinatoric_ratelaw = true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated stochastic chemical kinetics model SSAs for the reaction. Note,
for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X^2*Y*Z`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * binomial(X,2) *
binomial(Y,3)`.

Notes:
- Allocates
- `combinatoric_ratelaw = true` uses binomials in calculating the rate law, i.e. for `2S ->
 0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaw = false` then
 the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling
 factor.
"""
function jumpratelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate

    # if the `only_use_rate` argument is unused, update the rate according to the stoichiometries
    if !only_use_rate
        # Initiates the coefficient as "1", using the same type as in `substoich`.
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1

        # for each substrate, update the rate law (and coefficient) terms according to the
        # substrate's stoichiometry
        for (i, stoich) in enumerate(substoich)
            spec = substrates[i]

            # For symbolic substrates (e.g. parameters), simply update the rate law with the
            # appropriate expression. If not, we can simplify the expression (and also treat
            # the case with stoichiometry = 1 separately).
            if stoich isa Symbolic
                rl *= combinatoric_ratelaw ? binomial(spec, stoich) :
                      factorial(spec) / factorial(spec - stoich)
            else
                rl *= spec
                # Since `stoich` is a number, (if it is not 1) we add the binomial's denominator
                # to the `coef` (which is a number) directly. For Symbolic stoichiometry everything
                # is instead added to the rate law.
                isone(stoich) && continue
                for i in one(stoich):(stoich - one(stoich))
                    rl *= (spec - i)
                end
                combinatoric_ratelaw && (coef *= factorial(stoich))
            end
        end

        # if we use a combinatoric rate law, we have generated a numeric coefficient which we
        # will (if it is not 1) divide our rate law with (to generate the full rate law)
        combinatoric_ratelaw && !isequal(coef, one(coef)) && (rl /= coef)
    end
    rl
end

# if haveivdep = false then time-dependent rates will still be classified as mass action
"""
```julia
ismassaction(rx, rs; rxvars = get_variables(rx.rate),
 haveivdep = nothing,
 unknownset = Set(unknowns(rs)),
 ivset = nothing)
```

True if a given reaction is of mass action form, i.e. `rx.rate` does not depend
on any chemical species that correspond to unknowns of the system, and does not
depend explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as
 possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field
 explicitly depends on any independent variable (i.e. t or for spatial systems x,y,etc).
 If not set, will be automatically calculated.
- Optional: `unknownset`, set of unknowns which if the rxvars are within mean rx is
 non-mass action.
- Optional: `ivset`, a `Set` of the independent variables of the system. If not provided and
 the system is spatial, i.e. `isspatial(rs) == true`, it will be created with all the
 spatial variables and the time variable. If the rate expression contains any element of
 `ivset`, then `ismassaction(rx,rs) == false`. Pass a custom set to control this behavior.

Notes:
- Non-integer stoichiometry is treated as non-mass action. This includes
 symbolic variables/terms or floating point numbers for stoichiometric
 coefficients.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
    haveivdep::Union{Nothing, Bool} = nothing,
    unknownset = Set(get_unknowns(rs)), ivset = nothing)

    # we define non-integer (i.e. float or symbolic) stoich to be non-mass action
    ((eltype(rx.substoich) <: Integer) && (eltype(rx.prodstoich) <: Integer)) ||
        return false

    # if no dependencies must be zero order
    (length(rxvars) == 0) && return true

    # Reactions that depend on independent variables are not mass action. If iv dependency
    # was not marked when `ismassaction` was called, we check if the reaction depend on any ivs.
    if (haveivdep === nothing)
        if isspatial(rs)
            if (ivset === nothing)
                ivs = Set(get_sivs(rs))
                push!(ivs, get_iv(rs))
                ivdep = any(var -> var ∈ ivs, rxvars)
            else
                ivdep = any(var -> var ∈ ivset, rxvars)
            end
        else
            ivdep = any(var -> isequal(get_iv(rs), var), rxvars)
        end
            ivdep && return false
    else
        haveivdep && return false
    end

    # the `only_use_rate` option disables mass action for the reaction (e.g. if the  => arrow was
    # used in the DSL)
    rx.only_use_rate && return false

    # not mass action if have a non-constant, variable species in the rate expression
    @inbounds for var in rxvars
        (var in unknownset) && return false
    end

    return true
end

# from a reaction, create a `MassActionJump` structure (as used by JumpPrcoesses)
@inline function makemajump(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, netstoich = rx

    # if nothing changes in the reaction (e.g. `X --> X`), throws an error
    net_stoch = filter(p -> !drop_dynamics(p[1]), netstoich)
    isempty(net_stoch) &&
        error("$rx has no net stoichiometry change once accounting for constant and boundary condition species. This is not supported.")

    # creates the reactant stoichiometry vector (which maps each substrate to its stoichiometry)
    reactant_stoch = Vector{Pair{Any, eltype(substoich)}}()
    @inbounds for (i, spec) in enumerate(substrates)
        # move constant species into the rate
        if isconstant(spec)
            rate *= spec
            isone(substoich[i]) && continue
            for i in 1:(substoich[i] - 1)
                rate *= spec - i
            end
        else
            push!(reactant_stoch, substrates[i] => substoich[i])
        end
    end

    # for non-zeroth order reaction (e.g. `X1 --> X2`), and if we use a combinatorial rate law,
    # modifies the rate with its mass action coefficient (no need if coefficient is 1)
    zeroorder = (length(substoich) == 0)
    if (!zeroorder) && combinatoric_ratelaw
        coef = prod(factorial, substoich)
        (!isone(coef)) && (rate /= coef)
    end

    MassActionJump(Num(rate), reactant_stoch, net_stoch, scale_rates = false,
        useiszero = false)
end

# recursively visit each neighbor's rooted tree and mark everything in it as vrj
function dfs_mark!(isvrjvec, visited, depgraph, i)
    # Marks the current node as visited. Finds all neighbours.
    visited[i] = true
    nhbrs = depgraph[i]

    # for each neighbour to the input node, if we have not visited it, mark it as visited and 
    # recursively calls `dfs_mark!` on the neighbour.
    for nhbr in nhbrs
        if !visited[nhbr]
            isvrjvec[nhbr] = true
            dfs_mark!(isvrjvec, visited, depgraph, nhbr)
        end
    end
    nothing
end

# get_depgraph(rs)[i] is the list of reactions with rates depending on species changed by
# i'th reaction
function get_depgraph(rs)
    jdeps = asgraph(rs)
    vdeps = variable_dependencies(rs)
    eqeq_dependencies(jdeps, vdeps).fadjlist
end

# creates the full vector of jumps (which is then used to generate a JumpProcesses JumpSystem)
function assemble_jumps(rs; combinatoric_ratelaws = true)
    meqs = MassActionJump[]
    ceqs = ConstantRateJump[]
    veqs = VariableRateJump[]
    rxvars = []
    unknownset = Set(get_unknowns(rs))

    all(isspecies, unknownset) ||
        error("Conversion to JumpSystem currently requires all unknowns to be species.")
    isempty(get_rxs(rs)) &&
        error("Must give at least one reaction before constructing a JumpSystem.")

    # first we determine vrjs with an explicit time-dependent rate
    rxs = get_rxs(rs)
    isvrjvec = falses(length(rxs))
    havevrjs = false
    
    # for each reaction, check if any of its rate's variables is the time iv
    # this, and the next block, computes `isvrjvec` (which jumps are variable rate jumps) and 
    # `havevrjs` (whether the system has any variable rate jumps)
    for (i, rx) in enumerate(rxs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)
        @inbounds for rxvar in rxvars
            if isequal(rxvar, get_iv(rs))
                isvrjvec[i] = true
                havevrjs = true
                break
            end
        end
    end

    # now we determine vrj's that depend on species modified by a previous vrj
    if havevrjs
        depgraph = get_depgraph(rs)
        visited = falses(length(isvrjvec))
        for (i, isvrj) in enumerate(isvrjvec)
            if isvrj && !visited[i]
                # dfs from the vrj node to propagate vrj classification
                dfs_mark!(isvrjvec, visited, depgraph, i)
            end
        end
    end

    # loops through all reactions.
    # determines its type, creates the corresponding jump, and adds it to the correct vector.
    for (i, rx) in enumerate(rxs)
        # empties the `rxvars` cache and fills it with the variables in the reaction's rate
        # does it in this way to avoid allocation (to make things run faster)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)

        isvrj = isvrjvec[i]
        # if we are a mass action jump, create a mass action jump using `makemajump`
        if (!isvrj) && ismassaction(rx, rs; rxvars = rxvars, haveivdep = false,
                        unknownset = unknownset)
            push!(meqs, makemajump(rx; combinatoric_ratelaw = combinatoric_ratelaws))
        else
            # computes the jump rate and jump affect functions (these are identical for vrjs and crjs)
            rl = jumpratelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
            affect = Vector{Equation}()
            for (spec, stoich) in rx.netstoich
                # don't change species that are constant or BCs
                (!drop_dynamics(spec)) && push!(affect, spec ~ spec + stoich)
            end

            if isvrj
                push!(veqs, VariableRateJump(rl, affect))
            else
                push!(ceqs, ConstantRateJump(rl, affect))
            end
        end
    end
    vcat(meqs, ceqs, veqs)
end

### Equation Coupling ###

# adds any normal (differential or algebraic) equations to the equations (`eqs`) generated
# also adds any boundary value species to the unknowns vector
# also removes any equations produced by constant rate species
# also handles the addition of equations from boundary condition species
# also handles removing BC and constant species
# also handles the elimination of conserved quantities (if this was requested)
function add_coupled_eqs!(eqs, rs::ReactionSystem, ind_us, ind_sps; remove_conserved = false)
    # if there are BC species, put them after the independent species
    rs_us = get_unknowns(rs)
    us = any(isbc, rs_us) ? vcat(ind_us, filter(isbc, rs_us)) : ind_us
    ps = get_ps(rs)

    # if conservation law removal is designated, adds conservation law constants as parameters
    # and creates observables of the elimianted species.
    if remove_conserved
        nps = get_networkproperties(rs)
        
        # make dependent species observables and add conservation constants as parameters
        ps = vcat(ps, collect(eq.lhs for eq in nps.constantdefs))
        defs = copy(MT.defaults(rs))
        for eq in nps.constantdefs
            defs[eq.lhs] = eq.rhs
        end

        # adds the (by conservation laws) eliminated species as observables
        obs = copy(MT.observed(rs))
        append!(obs, nps.conservedeqs)
    else
        # if we have not eliminated any conservation laws, the defaults and observables are unchanged
        defs = MT.defaults(rs)
        obs = MT.observed(rs)
    end

    # Appends any coupled (differential and algebraic) equations from the systems to the
    # equations vector. Gives a note if coupled equations are combined with conservation law elimination.
    ceqs = Equation[eq for eq in nonreactions(rs)]
    if !isempty(ceqs)
        if remove_conserved
            @info """
                  Be careful mixing coupled equations and elimination of conservation laws.
                  Catalyst does not check that the conserved equations still hold for the
                  final coupled system of equations. Consider using `remove_conserved =
                  false` and instead calling ModelingToolkit.structural_simplify to simplify
                  any generated ODESystem or NonlinearSystem.
                  """
        end
        append!(eqs, ceqs)
    end

    eqs, us, ps, obs, defs
end

# used by flattened systems that don't support coupled equations currently
function error_if_coupled(::Type{T}, sys::ReactionSystem) where {T <: MT.AbstractSystem}
    any(eq -> eq isa Equation, get_eqs(sys)) &&
        error("Can not convert to a system of type ", T,
              " when there are coupled equations.")
    nothing
end

### Utility ###

# Throws an error when attempting to convert a spatial system to an unsuported type.
function spatial_convert_err(rs::ReactionSystem, systype)
    isspatial(rs) && error("Conversion to $systype is not supported for spatial networks.")
end

# Finds and differentials in an expression, and sets these to 0.
function remove_diffs(expr)
    if hasnode(Symbolics.is_derivative, expr)
        return replacenode(expr, diff_2_zero)
    else
        return expr
    end
end
diff_2_zero(expr) = (Symbolics.is_derivative(expr) ? 0 : expr)

COMPLETENESS_ERROR = "A ReactionSystem must be complete before it can be converted to other system types. A ReactionSystem can be marked as complete using the `complete` function."

# Used to, when required, display a warning about conservation law removal and remake.
function check_cons_warning(remove_conserved, remove_conserved_warn)
    (remove_conserved && remove_conserved_warn) || return
    @warn "You are creating a system or problem while eliminating conserved quantities. Please note,
        due to limitations / design choices in ModelingToolkit if you use the created system to 
        create a problem (e.g. an `ODEProblem`), or are directly creating a problem, you *should not*
        modify that problem's initial conditions for species (e.g. using `remake`). Changing initial 
        conditions must be done by creating a new Problem from your reaction system or the 
        ModelingToolkit system you converted it into with the new initial condition map. 
        Modification of parameter values is still possible, *except* for the modification of any
        conservation law constants ($CONSERVED_CONSTANT_SYMBOL), which is not possible. You might 
        get this warning when creating a problem directly.

        You can remove this warning by setting `remove_conserved_warn = false`."
end

### System Conversions ###

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.ODESystem`.

Keyword args and default values:
- `combinatoric_ratelaws = true` uses factorial scaling factors in calculating the rate law,
 i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
 `combinatoric_ratelaws = false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
 ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
 itself defaults to true).
- `remove_conserved = false`, if set to `true` will calculate conservation laws of the
 underlying set of reactions (ignoring coupled equations), and then apply them to reduce
 the number of equations.
"""
function Base.convert(::Type{<:ODESystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, remove_conserved = false, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, ODESystem)

    # Flattens the system. Next (if `remove_conserved = true`) uses `conservationlaws` to compute
    # these (these are then automatically cached in `rs`).
    flatrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(flatrs)

    # Find independent unknowns and species. Assembles the drift ODE, and then the full set of equations.
    # Here, `add_coupled_eqs!` performs some assembly in addition to adding coupled equations.
    ind_us, ind_sps = get_indep_us(flatrs, remove_conserved)
    eqs = assemble_drift(flatrs, ind_sps; combinatoric_ratelaws, remove_conserved,
                         include_zero_odes)
    eqs, us, ps, obs, defs = add_coupled_eqs!(eqs, flatrs, ind_us, ind_sps; remove_conserved)

    ODESystem(eqs, get_iv(flatrs), us, ps;
              observed = obs,
              name,
              defaults = _merge(defaults, defs),
              checks,
              continuous_events = MT.get_continuous_events(flatrs),
              discrete_events = MT.get_discrete_events(flatrs),
              kwargs...)
end

"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.NonlinearSystem`.

Keyword args and default values:
- `combinatoric_ratelaws = true` uses factorial scaling factors in calculating the rate law,
 i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
 `combinatoric_ratelaws = false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
 ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
 itself defaults to true).
- `remove_conserved = false`, if set to `true` will calculate conservation laws of the
 underlying set of reactions (ignoring coupled equations), and then apply them to reduce
 the number of equations.
"""
function Base.convert(::Type{<:NonlinearSystem}, rs::ReactionSystem; name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, remove_conserved = false, checks = false,
        remove_conserved_warn = true, default_u0 = Dict(), default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        all_differentials_permitted = false, kwargs...)
    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, NonlinearSystem)
    check_cons_warning(remove_conserved, remove_conserved_warn)
    if !isautonomous(rs)
        error("Attempting to convert a non-autonomous `ReactionSystem` (e.g. where some rate depend on $(get_iv(rs))) to a `NonlinearSystem`. This is not possible. if you are intending to compute system steady states, consider creating and solving a `SteadyStateProblem.")
    end

    # Flattens the system. Next (if `remove_conserved = true`) uses `conservationlaws` to compute
    # these (and automatically cache these in `rs`).
    flatrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(flatrs)

    # Find independent unknowns and species. Assembles the drift ODE, and then the full set of equations.
    # Here, `add_coupled_eqs!` performs some assembly in addition to adding coupled equations.
    ind_us, ind_sps = get_indep_us(flatrs, remove_conserved)
    eqs = assemble_drift(flatrs, ind_sps; combinatoric_ratelaws, remove_conserved,
                         as_odes = false, include_zero_odes)
    eqs, us, ps, obs, defs = add_coupled_eqs!(eqs, flatrs, ind_us, ind_sps; remove_conserved)

    # Throws a warning if there are differential equations in non-standard format.
    # Next, sets all differential terms to `0`.
    all_differentials_permitted || nonlinear_convert_differentials_check(rs)
    eqs = [remove_diffs(eq.lhs) ~ remove_diffs(eq.rhs) for eq in eqs]

    NonlinearSystem(eqs, us, ps;
                    name,
                    observed = obs,
                    defaults = _merge(defaults, defs),
                    checks,
                    kwargs...)
end

# Ideally, when `ReactionSystem`s are converted to `NonlinearSystem`s, any coupled ODEs should be 
# on the form D(X) ~ ..., where lhs is the time derivative w.r.t. a single variable, and the rhs
# does not contain any differentials. If this is not the case, we throw a warning to let the user
# know that they should be careful.
function nonlinear_convert_differentials_check(rs::ReactionSystem)
    for eq in filter(is_diff_equation, equations(rs))
        # For each differential equation, checks in order: 
        # If there is a differential on the right-hand side.
        # If the lhs is not on the form D(...).
        # If the lhs upper-level function is not a differential w.r.t. time.
        # If the content of the differential is not a variable (and nothing more).
        # If either of this is a case, throw the warning.
        if Symbolics._occursin(Symbolics.is_derivative, eq.rhs) ||
                               !Symbolics.istree(eq.lhs) ||
                               !isequal(Symbolics.operation(eq.lhs), Differential(get_iv(rs))) || 
                               (length(arguments(eq.lhs)) != 1) ||
                               !any(isequal(arguments(eq.lhs)[1]), nonspecies(rs))
            error("You are attempting to convert a `ReactionSystem` coupled with differential equations to a `NonlinearSystem`. However, some of these differentials are not of the form `D(x) ~ ...` where: 
                  (1) The left-hand side is a differential of a single variable with respect to the time-independent variable, and
                  (2) The right-hand side does not contain any differentials.
                  This is generally not permitted.
                                  
                  If you still would like to perform this conversion, please use the `all_differentials_permitted = true` option. In this case, all differential will be set to `0`. 
                  However, it is recommended to proceed with caution to ensure that the produced nonlinear equation makes sense for your intended application."
            )
        end
    end
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.SDESystem`.

Notes:
- `combinatoric_ratelaws = true` uses factorial scaling factors in calculating the rate law,
 i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
 `combinatoric_ratelaws = false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
 ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
 itself defaults to true).
- `remove_conserved = false`, if set to `true` will calculate conservation laws of the
 underlying set of reactions (ignoring coupled equations), and then apply them to reduce
 the number of equations.
- Does not currently support `ReactionSystem`s that includes coupled algebraic or
 differential equations.
"""
function Base.convert(::Type{<:SDESystem}, rs::ReactionSystem;
                      name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, checks = false, remove_conserved = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, SDESystem)
    check_cons_warning(remove_conserved, remove_conserved_warn)

    # Flattens the system. Next (if `remove_conserved = true`) uses `conservationlaws` to compute
    # these (and automatically cache these in `rs`).
    flatrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(flatrs)

    # Find independent unknowns and species. Assembles the drift ODE, and then the full set of
    # equations. Also assemblies the drift (stochasticity) terms.
    # Here, `add_coupled_eqs!` performs some assembly in addition to adding coupled equations.
    ind_us, ind_sps = get_indep_us(flatrs, remove_conserved)
    eqs = assemble_drift(flatrs, ind_sps; combinatoric_ratelaws, include_zero_odes,
                         remove_conserved)
    noiseeqs = assemble_diffusion(flatrs, ind_us, ind_sps; combinatoric_ratelaws, remove_conserved)
    eqs, us, ps, obs, defs = add_coupled_eqs!(eqs, flatrs, ind_us, ind_sps; remove_conserved)

    SDESystem(eqs, noiseeqs, get_iv(flatrs), us, ps;
        observed = obs,
        name,
        defaults = defs,
        checks,
        continuous_events = MT.get_continuous_events(flatrs),
        discrete_events = MT.get_discrete_events(flatrs),
        kwargs...)
end

"""
```julia
Base.convert(::Type{<:JumpSystem}, rs::ReactionSystem; combinatoric_ratelaws = true)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.JumpSystem`.

Notes:
- `combinatoric_ratelaws = true` uses binomials in calculating the rate law, i.e. for `2S ->
 0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaws = false` then
 the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling factor.
 Defaults to the value given when the `ReactionSystem` was constructed (which itself
 defaults to true).
- Does not currently support `ReactionSystem`s that includes coupled algebraic or
 differential equations.
- Does not currently support continuous events as these are not supported by
 `ModelingToolkit.JumpSystems`.
"""
function Base.convert(::Type{<:JumpSystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      remove_conserved = nothing, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, JumpSystem)
    (remove_conserved !== nothing) &&
        throw(ArgumentError("Catalyst does not support removing conserved species when converting to JumpSystems."))

    # Flattens the system. Performs post-flattening errors and warning checks.
    flatrs = Catalyst.flatten(rs)
    error_if_coupled(JumpSystem, flatrs)
    (length(MT.continuous_events(flatrs)) > 0) &&
        (@warn "continuous_events will be dropped as they are not currently supported by JumpSystems.")

    # Computes species, parameters, and equations. Adds boundary condition species to species vector.
    _, ind_sps = get_indep_us(flatrs, remove_conserved)
    any(isbc, get_unknowns(flatrs)) && (ind_sps = vcat(ind_sps, filter(isbc, get_unknowns(flatrs))))
    ps = get_ps(flatrs)
    eqs = assemble_jumps(flatrs; combinatoric_ratelaws)
    
    JumpSystem(eqs, get_iv(flatrs), ind_sps, ps;
               observed = MT.observed(flatrs),
               name,
               defaults = _merge(defaults, MT.defaults(flatrs)),
               checks,
               discrete_events = MT.discrete_events(flatrs),
               kwargs...)
end

### Problems ###

# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan,
        p = DiffEqBase.NullParameters(), args...;
        check_length = false, name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, remove_conserved = false, remove_conserved_warn = true,
        checks = false, structural_simplify = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
        remove_conserved, remove_conserved_warn)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        osys = MT.structural_simplify(osys)
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify = true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end

    return ODEProblem(osys, u0map, tspan, pmap, args...; check_length, kwargs...)
end

# NonlinearProblem from AbstractReactionNetwork
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0,
        p = DiffEqBase.NullParameters(), args...;
        name = nameof(rs), include_zero_odes = true,
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, remove_conserved_warn = true, checks = false,
        check_length = false, all_differentials_permitted = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    nlsys = convert(NonlinearSystem, rs; name, combinatoric_ratelaws, include_zero_odes,
        checks, all_differentials_permitted, remove_conserved, remove_conserved_warn)
    nlsys = complete(nlsys)
    return NonlinearProblem(nlsys, u0map, pmap, args...; check_length,
        kwargs...)
end

# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                               include_zero_odes = true, checks = false, check_length = false,
                               remove_conserved = false, structural_simplify = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    sde_sys = convert(SDESystem, rs; name, combinatoric_ratelaws,
                      include_zero_odes, checks, remove_conserved)
    p_matrix = zeros(length(get_unknowns(sde_sys)), numreactions(rs))

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        sde_sys = MT.structural_simplify(sde_sys)
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify = true` within `ODEProblem` call.")
    else
        sde_sys = complete(sde_sys)
    end
    
    return SDEProblem(sde_sys, u0map, tspan, pmap, args...; check_length,
        noise_rate_prototype = p_matrix, kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0, tspan::Tuple,
        p = DiffEqBase.NullParameters(), args...;
        name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    jsys = complete(jsys)
    return DiscreteProblem(jsys, u0map, tspan, pmap, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function JumpProcesses.JumpProblem(rs::ReactionSystem, prob, aggregator, args...;
        name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        checks = false, kwargs...)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    jsys = complete(jsys)
    return JumpProblem(jsys, prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0,
        p = DiffEqBase.NullParameters(), args...;
        check_length = false, name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, remove_conserved_warn = true, include_zero_odes = true,
        checks = false, structural_simplify = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
        remove_conserved, remove_conserved_warn)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        osys = MT.structural_simplify(osys)
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify = true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end

    return SteadyStateProblem(osys, u0map, pmap, args...; check_length, kwargs...)
end

### Symbolic Variable/Symbol Conversions ###

# convert symbol of the form :sys.a.b.c to a symbolic a.b.c
function _symbol_to_var(sys, sym)
    # The input is a normal Symbol, which exists in the input system (subsystems not considered).
    if hasproperty(sys, sym)
        var = getproperty(sys, sym, namespace = false)
    else
        # If sys does not contain sym it is possible that sys.subsystem contains it (or a 
        # sub-sub system etc.). Here, "₊" splits the system and subsystem names in the symbol.
        # We thus split the symbol to find, and go down the subsystem tree to attempt to find
        # the symbol. 
        strs = split(String(sym), "₊")   # need to check if this should be split of not!!!
        if length(strs) > 1
            var = getproperty(sys, Symbol(strs[1]), namespace = false)
            for str in view(strs, 2:length(strs))
                var = getproperty(var, Symbol(str), namespace = true)
            end
        else
            throw(ArgumentError("System $(nameof(sys)): variable $sym does not exist"))
        end
    end
    var
end

"""
 symmap_to_varmap(sys, symmap)

Given a system and map of `Symbol`s to values, generates a map from
corresponding symbolic variables/parameters to the values that can be used to
pass initial conditions and parameter mappings.

For example,
```julia
sir = @reaction_network sir begin
 β, S + I --> 2I
 ν, I --> R
end
subsys = @reaction_network subsys begin
 k, A --> B
end
@named sys = compose(sir, [subsys])
```
gives
```
Model sys with 3 equations
Unknowns (5):
 S(t)
 I(t)
 R(t)
 subsys₊A(t)
 subsys₊B(t)
Parameters (3):
 β
 ν
 subsys₊k
```
to specify initial condition and parameter mappings from *symbols* we can use
```julia
symmap = [:S => 1.0, :I => 1.0, :R => 1.0, :subsys₊A => 1.0, :subsys₊B => 1.0]
u0map  = symmap_to_varmap(sys, symmap)
pmap   = symmap_to_varmap(sys, [:β => 1.0, :ν => 1.0, :subsys₊k => 1.0])
```
`u0map` and `pmap` can then be used as input to various problem types.

Notes:
- Any `Symbol`, `sym`, within `symmap` must be a valid field of `sys`. i.e.
 `sys.sym` must be defined.
"""
function symmap_to_varmap(sys, symmap::Tuple)
    # If all inputs are of the form `:sym => value`, apply `_symbol_to_var` function. If not,
    # pass input through unchanged.
    if all(p -> p isa Pair{Symbol}, symmap)
        return ((_symbol_to_var(sys, sym) => val for (sym, val) in symmap)...,)
    else
        return symmap
    end
end

function symmap_to_varmap(sys, symmap::AbstractArray{Pair{Symbol, T}}) where {T}
    [_symbol_to_var(sys, sym) => val for (sym, val) in symmap]
end

function symmap_to_varmap(sys, symmap::Dict{Symbol, T}) where {T}
    Dict(_symbol_to_var(sys, sym) => val for (sym, val) in symmap)
end

# don't permute any other types and let varmap_to_vars handle erroring
symmap_to_varmap(sys, symmap) = symmap
#error("symmap_to_varmap requires a Dict, AbstractArray or Tuple to map Symbols to values.")

### Other Conversion-related Functions ###

# the following function is adapted from SymbolicUtils.jl v.19
# later on (Spetember 2023) modified by Torkel and Shashi (now assumes input not on polynomial form, which is handled elsewhere, the previous version failed in these cases anyway).
# Copyright (c) 2020: Shashi Gowda, Yingbo Ma, Mason Protter, Julia Computing.
# MIT license
"""
 to_multivariate_poly(polyeqs::AbstractVector{BasicSymbolic{Real}})

Convert the given system of polynomial equations to multivariate polynomial representation.
For example, this can be used in HomotopyContinuation.jl functions.
"""
function to_multivariate_poly(polyeqs::AbstractVector{Symbolics.BasicSymbolic{Real}})
    @assert length(polyeqs) >= 1 "At least one expression must be passed to `multivariate_poly`."

    pvar2sym, sym2term = SymbolicUtils.get_pvar2sym(), SymbolicUtils.get_sym2term()
    ps = map(polyeqs) do x
        if iscall(x) && operation(x) == (/)
            error("We should not be able to get here, please contact the package authors.")
        else
            PolyForm(x, pvar2sym, sym2term).p
        end
    end

    ps
end
