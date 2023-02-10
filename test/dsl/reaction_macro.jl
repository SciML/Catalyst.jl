#! format: off

# Fetch pakages.
using Catalyst, ModelingToolkit

# Tests that weird arrows work.

@test isequal((@reaction k, 0 --> X), (@reaction k, X <-- 0))
@test isequal((@reaction k, 0 --> X), (@reaction k, X ⟻ 0))
@test isequal((@reaction k, 0 --> X), (@reaction k, 0 → X))
@test isequal((@reaction k, 0 --> X), (@reaction k, 0 ⥟ X))

# Tests disabling law of mass action.
@test isequal((@reaction k, X --> 0), (@reaction k*X, X => 0))