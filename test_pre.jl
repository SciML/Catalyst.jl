using Catalyst
using ModelingToolkit: Pre

# Set up variables
t = default_t()

# Test if Pre() operator works with discrete events
@species X(t)
@parameters p

# Test basic discrete event with Pre() operator
try
    rn = @reaction_network begin
        @discrete_events [1.0] => [X ~ Pre(X) + 1]
        p, 0 --> X
    end
    println("✓ Pre() operator works with discrete events")
    
    # Test continuous event with Pre() operator  
    rn2 = @reaction_network begin
        @continuous_events [X ~ 2.0] => [X ~ Pre(X) - 0.5]
        p, 0 --> X
    end
    println("✓ Pre() operator works with continuous events")

catch e
    println("✗ Error with Pre() operator: ", e)
end