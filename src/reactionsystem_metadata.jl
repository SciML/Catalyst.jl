### System-Level Metadata Key Types and Accessors ###

# This file defines metadata key types that can be stored in a ReactionSystem's
# `metadata` field (a `Base.ImmutableDict{DataType, Any}`), along with convenience
# accessor functions. These are primarily intended for use by file parsers that
# need to attach extra mapping data to a system.

"""
    U0Map

Metadata key for storing initial condition / species value mappings on a
[`ReactionSystem`](@ref). Intended for use by file parsers that need to preserve
mappings in a format different from `initial_conditions`.

See also: [`has_u0_map`](@ref), [`get_u0_map`](@ref), [`set_u0_map`](@ref)
"""
struct U0Map end

"""
    ParameterMap

Metadata key for storing parameter value mappings on a [`ReactionSystem`](@ref).
Intended for use by file parsers that need to preserve parameter mappings in a
format different from `initial_conditions`.

See also: [`has_parameter_map`](@ref), [`get_parameter_map`](@ref), [`set_parameter_map`](@ref)
"""
struct ParameterMap end

"""
    has_u0_map(rs::ReactionSystem)

Returns `true` if the [`ReactionSystem`](@ref) has a `U0Map` metadata entry.
"""
has_u0_map(rs::ReactionSystem) = hasmetadata(rs, U0Map)

"""
    get_u0_map(rs::ReactionSystem)

Returns the `U0Map` metadata from the [`ReactionSystem`](@ref), or `nothing` if
no `U0Map` has been set.
"""
get_u0_map(rs::ReactionSystem) = getmetadata(rs, U0Map, nothing)

"""
    set_u0_map(rs::ReactionSystem, u0map)

Returns a **new** [`ReactionSystem`](@ref) with the `U0Map` metadata set to `u0map`.
The original system is not modified.
"""
set_u0_map(rs::ReactionSystem, u0map) = setmetadata(rs, U0Map, u0map)

"""
    has_parameter_map(rs::ReactionSystem)

Returns `true` if the [`ReactionSystem`](@ref) has a `ParameterMap` metadata entry.
"""
has_parameter_map(rs::ReactionSystem) = hasmetadata(rs, ParameterMap)

"""
    get_parameter_map(rs::ReactionSystem)

Returns the `ParameterMap` metadata from the [`ReactionSystem`](@ref), or `nothing` if
no `ParameterMap` has been set.
"""
get_parameter_map(rs::ReactionSystem) = getmetadata(rs, ParameterMap, nothing)

"""
    set_parameter_map(rs::ReactionSystem, pmap)

Returns a **new** [`ReactionSystem`](@ref) with the `ParameterMap` metadata set to `pmap`.
The original system is not modified.
"""
set_parameter_map(rs::ReactionSystem, pmap) = setmetadata(rs, ParameterMap, pmap)
