using Catalyst, Test

rn = @reaction_network begin
    (1, 2), A + B <--> C 
    (3, 2), D <--> E 
    (0.1, 0.2), E <--> F
    6, F --> G
    7, H --> G
    5, 0 --> K
    3, A + B --> Z + C
end

S = netstoichmat(rn)
C = conservationlaws(S)

@test size(C,1) == 3

b = [ 0, 0, 0, 1, 1, 1, 1, 1, 0, 0 ]
@test any(b == C[i,:] for i in 1:size(C,1))

# For the A + B <--> C subsystem one of these must occur
# as a conservation law
D = [ 1 -1 0 0 0 0 0 0 0 0;
      -1 1 0 0 0 0 0 0 0 0
      0  1 1 0 0 0 0 0 0 0 ]

@test any(D[j,:] == C[i,:] for i in 1:size(C,1), j in size(D,1))
