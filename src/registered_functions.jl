using Symbolics

# Registers the Michaelis-Menten function.
mm(X,v,K) = v*X / (X + K)

@register mm(X,v,K);

Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{1}) = (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{2}) = args[1]/(args[1]+args[3])
Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{3}) = - args[2]*args[1]/(args[1]+args[3])^2

# Registers the repressing Michaelis-Menten function.
mmr(X,v,K) = v*K / (X + K)

@register mmr(X,v,K); 

Symbolics.derivative(::typeof(mmr), args::NTuple{3,Any}, ::Val{1}) = - (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(mmr), args::NTuple{3,Any}, ::Val{2}) = args[3]/(args[1]+args[3])
Symbolics.derivative(::typeof(mmr), args::NTuple{3,Any}, ::Val{3}) = args[2]*args[1]/(args[1]+args[3])^2


# Registers the Hill function.
hill(X,v,K,n) = v*(X^n) / (X^n + K^n)

@register hill(X,v,K,n);
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{1}) = args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{2}) = (args[1]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{3}) = - args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[1])-log(args[3]))  /  (args[1]^args[4] + args[3]^args[4])^2

# Registers the repressing hill function (alterantive to negative n).
hillr(X,v,K,n) = v*(K^n) / (X^n + K^n)

@register hillr(X,v,K,n);

Symbolics.derivative(::typeof(hillr), args::NTuple{4,Any}, ::Val{1}) = - args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hillr), args::NTuple{4,Any}, ::Val{2}) = (args[3]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(hillr), args::NTuple{4,Any}, ::Val{3}) = args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hillr), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[3])-log(args[1]))  /  (args[1]^args[4] + args[3]^args[4])^2

# Registers the activation/repressing hill function.
hillar(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)
@register hillar(X,Y,v,K,n);

Symbolics.derivative(::typeof(hillar), args::NTuple{5,Any}, ::Val{1}) = args[3] * args[5] * (args[1]^(args[5]-1)) * (args[2]^args[5]+args[4]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillar), args::NTuple{5,Any}, ::Val{2}) = - args[3] * args[5] * (args[2]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillar), args::NTuple{5,Any}, ::Val{3}) = (args[1]^args[5])   /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
Symbolics.derivative(::typeof(hillar), args::NTuple{5,Any}, ::Val{4}) = - args[3] * args[5] * (args[3]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillar), args::NTuple{5,Any}, ::Val{5}) = args[3] * (args[1]^args[5])  *  (log(args[1])*(args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5])*log(args[2]) - (args[4]^args[5])*log(args[4]))   /   (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
