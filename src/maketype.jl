function maketype(name, syms, scale_noise; params = Symbol[], 
                reactions = Vector{ReactionStruct}(undef, 0),
                syms_to_ints = OrderedDict{Symbol,Int}(),
                params_to_ints = OrderedDict{Symbol,Int}(),
                properties = Dict{Symbol,Any}())
    typeex = :(mutable struct $name <: DiffEqBase.AbstractReactionNetwork
        syms::Vector{Symbol}
        scale_noise::Symbol
        params::Vector{Symbol}
        reactions::Vector{ReactionStruct}
        syms_to_ints::OrderedDict{Symbol,Int}
        params_to_ints::OrderedDict{Symbol,Int}
        properties::Dict{Symbol,Any}
    end)

    # Make the default constructor
    constructorex = :($(name)(;
                    $(Expr(:kw, :syms, syms)),
                    $(Expr(:kw, :scale_noise, Meta.quot(scale_noise))),
                    $(Expr(:kw, :params, params)),
                    $(Expr(:kw, :reactions, reactions)),
                    $(Expr(:kw, :syms_to_ints, syms_to_ints)),
                    $(Expr(:kw, :params_to_ints, params_to_ints)),
                    $(Expr(:kw, :properties, properties))) =
                    $(name)(syms, scale_noise, params, reactions, syms_to_ints, params_to_ints, properties)) |> esc

    # Make the type instance using the default constructor
    typeex, constructorex
end
