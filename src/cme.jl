
"""
Performs a breadth-first graph traversal of a  to enumerate all possible states.

The truncation kwarg should be given as a vector of Int, denoting the maximum species' counts.
If, upon taking traversing a reaction the amount of any species is greater than the value in the truncation vector, 
the graph traversal will not add the state to the list of possible states.
"""
function traverse_reactionsystem(rs, u0; truncation=nothing)
    num_molecules = sum(u0)
    num_species = length(u0)
    state_space = []
    Q = Queue{typeof(u0)}()
    enqueue!(Q, u0)
    push!(state_space, u0)
    nsm = netstoichmat(rs)
    edges = eachcol(nsm)
    while !isempty(Q)
        w = dequeue!(Q)
        for e in edges
            x = w .+ e
            if x ∉ state_space && all(y >= 0 for y in x)
                if truncation === nothing 
                    if sum(x) > num_molecules 
                        error("need to provide truncation for infinite reaction system")
                    else 
                        push!(state_space, x)
                        enqueue!(Q, x)
                    end
                else 
                    if all(x[i] <= truncation[i] for i in 1:num_species)
                        push!(state_space, x)
                        enqueue!(Q, x)
                    end
                end
            end
        end
    end
    return state_space
end

function ODESystem_cme(rs, u0; truncation=nothing, kwargs...)
	state_space = traverse_reactionsystem(rs, u0; truncation=truncation)

	N = length(state_space)
	d = Dict(state_space .=> 1:N)

    # initializes the state corresponding to u0 to have probability 1
	defs = zeros(N)
	start_state = d[u0]
	defs[start_state] = 1.


	@variables t u[1:N](t) = defs
	D = Differential(t)

	rxns = reactions(rs)
	M = numreactions(rs)
	jump_rates = jumpratelaw.(rxns)
	nsm = netstoichmat(rs)
	species = states(rs)
	
	eqs = Vector{Equation}(undef, N)
	for (i, x) in enumerate(state_space)
		eq = Num(0)
		for j in 1:M
			x′ = x .- @view(nsm[:, j]) # traverse a reaction and get the new counts
			p_x′ = haskey(d, x′) ? u[d[x′]] : 0 # if the new state is not in `state_space`, set probability to 0
			t1 = substitute(jump_rates[j], Dict(species .=> x′)) * p_x′
			t2 = substitute(jump_rates[j], Dict(species .=> x)) * u[d[x]]
			eq += t1 - t2
		end 
		eqs[i] = D(u[d[x]]) ~ eq
	end
	ODESystem(eqs, t, u, parameters(rs); name=nameof(rs), kwargs...)
end
