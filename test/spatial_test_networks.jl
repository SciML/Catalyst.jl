### Fetch packages ###
using Catalyst, Graphs

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

### Helper Functions ###

# Generates ranomised intiial condition or paraemter values.
rand_v_vals(grid) = rand(rng, nv(grid))
rand_v_vals(grid, x::Number) = rand_v_vals(grid) * x
rand_e_vals(grid) = rand(rng, ne(grid))
rand_e_vals(grid, x::Number) = rand_e_vals(grid) * x
function make_u0_matrix(value_map, vals, symbols)
    (length(symbols) == 0) && (return zeros(0, length(vals)))
    d = Dict(value_map)
    return [(d[s] isa Vector) ? d[s][v] : d[s] for s in symbols, v in 1:length(vals)]
end

# Gets a symbol list of spatial parameters.
function spatial_param_syms(lrs::LatticeReactionSystem)
    ModelingToolkit.getname.(diffusion_parameters(lrs))
end

# Converts to integer value (for JumpProcess simulations).
function make_values_int(values::Vector{<:Pair})
    [val[1] => round.(Int64, val[2]) for val in values]
end
make_values_int(values::Matrix{<:Number}) = round.(Int64, values)
make_values_int(values::Vector{<:Number}) = round.(Int64, values)
make_values_int(values::Vector{Vector}) = [round.(Int64, vals) for vals in values]

### Declares Models ###

# Small non-stiff system.
SIR_system = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
SIR_p = [:α => 0.1 / 1000, :β => 0.01]
SIR_u0 = [:S => 999.0, :I => 1.0, :R => 0.0]

SIR_tr_S = @transport_reaction dS S
SIR_tr_I = @transport_reaction dI I
SIR_tr_R = @transport_reaction dR R
SIR_srs_1 = [SIR_tr_S]
SIR_srs_2 = [SIR_tr_S, SIR_tr_I, SIR_tr_R]

# Small non-stiff system.
binding_system = @reaction_network begin (k1, k2), X + Y <--> XY end
binding_tr_X = @transport_reaction dX X
binding_tr_Y = @transport_reaction dY Y
binding_tr_XY = @transport_reaction dXY XY
binding_srs = [binding_tr_X, binding_tr_Y, binding_tr_XY]
binding_u0 = [:X => 1.0, :Y => 2.0, :XY => 0.5]
binding_p = [:k1 => 2.0, :k2 => 0.1, :dX => 3.0, :dY => 5.0, :dXY => 2.0]

# Mid-sized non-stiff system.
CuH_Amination_system = @reaction_network begin
    10.0^kp1, CuoAc + Ligand --> CuoAcLigand
    10.0^kp2, CuoAcLigand + Silane --> CuHLigand + SilaneOAc
    10.0^k1, CuHLigand + Styrene --> AlkylCuLigand
    10.0^k_1, AlkylCuLigand --> CuHLigand + Styrene
    10.0^k2, AlkylCuLigand + Amine_E --> AlkylAmine + Cu_ELigand
    10.0^k_2, AlkylAmine + Cu_ELigand --> AlkylCuLigand + Amine_E
    10.0^k3, Cu_ELigand + Silane --> CuHLigand + E_Silane
    10.0^kam, CuHLigand + Amine_E --> Amine + Cu_ELigand
    10.0^kdc, CuHLigand + CuHLigand --> Decomposition
end
CuH_Amination_p = [
    :kp1 => 1.2,
    :kp2 => -0.72,
    :k1 => 0.57,
    :k_1 => -3.5,
    :k2 => -0.35,
    :k_2 => -0.77,
    :k3 => -0.025,
    :kam => -2.6,
    :kdc => -3.0,
]
CuH_Amination_u0 = [
    :CuoAc => 0.0065,
    :Ligand => 0.0072,
    :CuoAcLigand => 0.0,
    :Silane => 0.65,
    :CuHLigand => 0.0,
    :SilaneOAc => 0.0,
    :Styrene => 0.16,
    :AlkylCuLigand => 0.0,
    :Amine_E => 0.39,
    :AlkylAmine => 0.0,
    :Cu_ELigand => 0.0,
    :E_Silane => 0.0,
    :Amine => 0.0,
    :Decomposition => 0.0,
]

CuH_Amination_tr_1 = @transport_reaction D1 CuoAc
CuH_Amination_tr_2 = @transport_reaction D2 Silane
CuH_Amination_tr_3 = @transport_reaction D3 Cu_ELigand
CuH_Amination_tr_4 = @transport_reaction D4 Amine
CuH_Amination_tr_5 = @transport_reaction D5 CuHLigand
CuH_Amination_srs_1 = [CuH_Amination_tr_1]
CuH_Amination_srs_2 = [
    CuH_Amination_tr_1,
    CuH_Amination_tr_2,
    CuH_Amination_tr_3,
    CuH_Amination_tr_4,
    CuH_Amination_tr_5,
]

# Small stiff system.
brusselator_system = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end
brusselator_p = [:A => 1.0, :B => 4.0]

brusselator_tr_x = @transport_reaction dX X
brusselator_tr_y = @transport_reaction dY Y
brusselator_srs_1 = [brusselator_tr_x]
brusselator_srs_2 = [brusselator_tr_x, brusselator_tr_y]

# Mid-sized stiff system.
# Unsure about stifness, but non-spatial version oscillates for this parameter set.
sigmaB_system = @reaction_network begin
    kDeg, (w, w2, w2v, v, w2v2, vP, σB, w2σB) ⟶ ∅
    kDeg, vPp ⟶ phos
    (kBw, kDw), 2w ⟷ w2
    (kB1, kD1), w2 + v ⟷ w2v
    (kB2, kD2), w2v + v ⟷ w2v2
    kK1, w2v ⟶ w2 + vP
    kK2, w2v2 ⟶ w2v + vP
    (kB3, kD3), w2 + σB ⟷ w2σB
    (kB4, kD4), w2σB + v ⟷ w2v + σB
    (kB5, kD5), vP + phos ⟷ vPp
    kP, vPp ⟶ v + phos
    v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ σB
    λW * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ w
    λV * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ v
end
sigmaB_p = [:kBw => 3600, :kDw => 18, :kB1 => 3600, :kB2 => 3600, :kB3 => 3600,
    :kB4 => 1800, :kB5 => 3600,
    :kD1 => 18, :kD2 => 18, :kD3 => 18, :kD4 => 1800, :kD5 => 18, :kK1 => 36, :kK2 => 12,
    :kP => 180, :kDeg => 0.7,
    :v0 => 0.4, :F => 30, :K => 0.2, :λW => 4, :λV => 4.5]
sigmaB_u0 = [
    :w => 1.0,
    :w2 => 1.0,
    :w2v => 1.0,
    :v => 1.0,
    :w2v2 => 1.0,
    :vP => 1.0,
    :σB => 1.0,
    :w2σB => 1.0,
    :vPp => 0.0,
    :phos => 0.4,
]

sigmaB_tr_σB = @transport_reaction DσB σB
sigmaB_tr_w = @transport_reaction Dw w
sigmaB_tr_v = @transport_reaction Dv v
sigmaB_srs_1 = [sigmaB_tr_σB]
sigmaB_srs_2 = [sigmaB_tr_σB, sigmaB_tr_w, sigmaB_tr_v]

### Declares Lattices ###

# Grids.
very_small_2d_grid = Graphs.grid([2, 2])
small_2d_grid = Graphs.grid([5, 5])
medium_2d_grid = Graphs.grid([20, 20])
large_2d_grid = Graphs.grid([100, 100])

small_3d_grid = Graphs.grid([5, 5, 5])
medium_3d_grid = Graphs.grid([20, 20, 20])
large_3d_grid = Graphs.grid([100, 100, 100])

# Paths.
short_path = path_graph(100)
long_path = path_graph(1000)

# Unconnected graphs.
unconnected_graph = SimpleGraph(10)

# Undirected cycle.
undirected_cycle = cycle_graph(49)

# Directed cycle.
small_directed_cycle = cycle_graph(100)
large_directed_cycle = cycle_graph(1000)