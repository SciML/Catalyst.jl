using Catalyst, Test

rn = @reaction_network arrowtest begin
    @parameters a1 a2 k1 k2 b1
    (a1, a2), C <--> 0
    (k1, k2), A + B <--> C
    b1, 0 <-- B
end

rn2 = @reaction_network arrowtest begin
    @parameters a1 a2 k1 k2 b1
    a1, C --> 0
    a2, 0 --> C
    k1, A + B --> C
    k2, C --> A + B
    b1, B --> 0
end

@test rn == rn2
