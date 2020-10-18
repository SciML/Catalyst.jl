# The Reaction DSL - Advanced
This section covers some of the more advanced syntax and features for building
chemical reaction network models (still not very complicated!).

#### User defined functions in reaction rates
The reaction network DSL cannot "see" user defined functions. E.g., this is not
correct syntax:

```julia
myHill(x) = 2.0*x^3/(x^3+1.5^3)
rn = @reaction_network begin
  myHill(X), ∅ → X
end
```
However, it is possible to define functions in such a way that the DSL can see
them using the `@reaction_func` macro:

```julia
@reaction_func myHill(x) = 2.0*x^3/(x^3+1.5^3)
rn = @reaction_network begin
  myHill(X), ∅ → X
end
```

#### Ignoring mass action kinetics
While generally one wants the reaction rate to use the law of mass action, so
the reaction
```julia
rn = @reaction_network begin
  k, X → ∅
end k
```
occurs at the rate ``d[X]/dt = -k[X]``, it is possible to ignore this by using
any of the following non-filled arrows when declaring the reaction: `⇐`, `⟽`,
`⇒`, `⟾`, `⇔`, `⟺`. This means that the reaction

```julia
rn = @reaction_network begin
  k, X ⇒ ∅
end k
```

will occur at rate ``d[X]/dt = -k`` (which might become a problem since ``[X]``
will be degraded at a constant rate even when very small or equal to 0).
