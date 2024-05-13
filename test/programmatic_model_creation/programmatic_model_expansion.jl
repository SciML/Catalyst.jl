#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, Test
using ModelingToolkit: get_ps, get_unknowns, get_eqs, get_systems, get_iv, getname, nameof

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Sets the default `t` to use.
t = default_t()

# Fetch test networks.
include("../test_networks.jl")

# Declares a helper function.
function unpacksys(sys)
    get_eqs(sys), get_iv(sys), get_ps(sys), nameof(sys), get_systems(sys)
end

### Basic Tests ###

# Tests construction of empty reaction networks.
let
    empty_network_1 = @reaction_network
    eqs, iv, ps, name, systems = unpacksys(empty_network_1)
    @test length(eqs) == 0
    @test nameof(iv) == :t
    @test length(get_unknowns(empty_network_1)) == 0
    @test length(ps) == 0
end
