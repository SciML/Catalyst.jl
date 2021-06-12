using Symbolics

# Registers the Michaelis-Menten function.
mm(X,v,K) = v*X / (X + K)
Mm(X,v,K) = v*X / (X + K)
MM(X,v,K) = v*X / (X + K)

@register mm(X,v,K); @register Mm(X,v,K); @register MM(X,v,K);

Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{1}) = (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{2}) = args[1]/(args[1]+args[3])
Symbolics.derivative(::typeof(mm), args::NTuple{3,Any}, ::Val{3}) = - args[2]*args[1]/(args[1]+args[3])^2

Symbolics.derivative(::typeof(Mm), args::NTuple{3,Any}, ::Val{1}) = (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(Mm), args::NTuple{3,Any}, ::Val{2}) = args[1]/(args[1]+args[3])
Symbolics.derivative(::typeof(Mm), args::NTuple{3,Any}, ::Val{3}) = - args[2]*args[1]/(args[1]+args[3])^2

Symbolics.derivative(::typeof(MM), args::NTuple{3,Any}, ::Val{1}) = (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(MM), args::NTuple{3,Any}, ::Val{2}) = args[1]/(args[1]+args[3])
Symbolics.derivative(::typeof(MM), args::NTuple{3,Any}, ::Val{3}) = - args[2]*args[1]/(args[1]+args[3])^2

# Registers the repressing Michaelis-Menten function.
mmR(X,v,K) = v*K / (X + K)
MmR(X,v,K) = v*K / (X + K)
MMR(X,v,K) = v*K / (X + K)

@register mmR(X,v,K); @register MmR(X,v,K); @register MMR(X,v,K);

Symbolics.derivative(::typeof(mmR), args::NTuple{3,Any}, ::Val{1}) = - (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(mmR), args::NTuple{3,Any}, ::Val{2}) = args[3]/(args[1]+args[3])
Symbolics.derivative(::typeof(mmR), args::NTuple{3,Any}, ::Val{3}) = args[2]*args[1]/(args[1]+args[3])^2

Symbolics.derivative(::typeof(MmR), args::NTuple{3,Any}, ::Val{1}) = - (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(MmR), args::NTuple{3,Any}, ::Val{2}) = args[3]/(args[1]+args[3])
Symbolics.derivative(::typeof(MmR), args::NTuple{3,Any}, ::Val{3}) = args[2]*args[1]/(args[1]+args[3])^2

Symbolics.derivative(::typeof(MMR), args::NTuple{3,Any}, ::Val{1}) = - (args[2]*args[3]) / (args[1]+args[3])^2
Symbolics.derivative(::typeof(MMR), args::NTuple{3,Any}, ::Val{2}) = args[3]/(args[1]+args[3])
Symbolics.derivative(::typeof(MMR), args::NTuple{3,Any}, ::Val{3}) = args[2]*args[1]/(args[1]+args[3])^2


# Registers the Hill function.
hill(X,v,K,n) = v*(X^n) / (X^n + K^n)
Hill(X,v,K,n) = v*(X^n) / (X^n + K^n)

@register hill(X,v,K,n); @register Hill(X,v,K,n);

Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{1}) = args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{2}) = (args[1]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{3}) = - args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hill), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[1])-log(args[3]))  /  (args[1]^args[4] + args[3]^args[4])^2

Symbolics.derivative(::typeof(Hill), args::NTuple{4,Any}, ::Val{1}) = args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(Hill), args::NTuple{4,Any}, ::Val{2}) = (args[1]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(Hill), args::NTuple{4,Any}, ::Val{3}) = - args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(Hill), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[1])-log(args[3]))  /  (args[1]^args[4] + args[3]^args[4])^2


# Registers the repressing hill function (alterantive to negative n).
hillR(X,v,K,n) = v*(K^n) / (X^n + K^n)
HillR(X,v,K,n) = v*(K^n) / (X^n + K^n)

@register hillR(X,v,K,n); @register HillR(X,v,K,n);

Symbolics.derivative(::typeof(hillR), args::NTuple{4,Any}, ::Val{1}) = - args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hillR), args::NTuple{4,Any}, ::Val{2}) = (args[3]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(hillR), args::NTuple{4,Any}, ::Val{3}) = args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(hillR), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[3])-log(args[1]))  /  (args[1]^args[4] + args[3]^args[4])^2

Symbolics.derivative(::typeof(HillR), args::NTuple{4,Any}, ::Val{1}) = - args[2] * args[4] * (args[3]^args[4]) * (args[1]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(HillR), args::NTuple{4,Any}, ::Val{2}) = (args[3]^args[4])  /  (args[1]^args[4] + args[3]^args[4])
Symbolics.derivative(::typeof(HillR), args::NTuple{4,Any}, ::Val{3}) = args[2] * args[4] * (args[1]^args[4]) * (args[3]^(args[4]-1))  /  (args[1]^args[4] + args[3]^args[4])^2
Symbolics.derivative(::typeof(HillR), args::NTuple{4,Any}, ::Val{4}) = args[2] * (args[1]^args[4]) * (args[3]^args[4]) * (log(args[3])-log(args[1]))  /  (args[1]^args[4] + args[3]^args[4])^2


# Registers the activation/repressing hill function.
hillC(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)
hillAR(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)
HillC(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)
HillAR(X,Y,v,K,n) = v*(X^n) / (X^n + Y^n + K^n)

@register hillC(X,Y,v,K,n); @register hillAR(X,Y,v,K,n); @register HillC(X,Y,v,K,n); @register HillAR(X,Y,v,K,n);

Symbolics.derivative(::typeof(hillC), args::NTuple{5,Any}, ::Val{1}) = args[3] * args[5] * (args[1]^(args[5]-1)) * (args[2]^args[5]+args[4]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillC), args::NTuple{5,Any}, ::Val{2}) = - args[3] * args[5] * (args[2]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillC), args::NTuple{5,Any}, ::Val{3}) = (args[1]^args[5])   /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
Symbolics.derivative(::typeof(hillC), args::NTuple{5,Any}, ::Val{4}) = - args[3] * args[5] * (args[3]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillC), args::NTuple{5,Any}, ::Val{5}) = args[3] * (args[1]^args[5])  *  (log(args[1])*(args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5])*log(args[2]) - (args[4]^args[5])*log(args[4]))   /   (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2

Symbolics.derivative(::typeof(hillAR), args::NTuple{5,Any}, ::Val{1}) = args[3] * args[5] * (args[1]^(args[5]-1)) * (args[2]^args[5]+args[4]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillAR), args::NTuple{5,Any}, ::Val{2}) = - args[3] * args[5] * (args[2]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillAR), args::NTuple{5,Any}, ::Val{3}) = (args[1]^args[5])   /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
Symbolics.derivative(::typeof(hillAR), args::NTuple{5,Any}, ::Val{4}) = - args[3] * args[5] * (args[3]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(hillAR), args::NTuple{5,Any}, ::Val{5}) = args[3] * (args[1]^args[5])  *  (log(args[1])*(args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5])*log(args[2]) - (args[4]^args[5])*log(args[4]))   /   (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2

Symbolics.derivative(::typeof(HillC), args::NTuple{5,Any}, ::Val{1}) = args[3] * args[5] * (args[1]^(args[5]-1)) * (args[2]^args[5]+args[4]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillC), args::NTuple{5,Any}, ::Val{2}) = - args[3] * args[5] * (args[2]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillC), args::NTuple{5,Any}, ::Val{3}) = (args[1]^args[5])   /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
Symbolics.derivative(::typeof(HillC), args::NTuple{5,Any}, ::Val{4}) = - args[3] * args[5] * (args[3]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillC), args::NTuple{5,Any}, ::Val{5}) = args[3] * (args[1]^args[5])  *  (log(args[1])*(args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5])*log(args[2]) - (args[4]^args[5])*log(args[4]))   /   (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2

Symbolics.derivative(::typeof(HillAR), args::NTuple{5,Any}, ::Val{1}) = args[3] * args[5] * (args[1]^(args[5]-1)) * (args[2]^args[5]+args[4]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillAR), args::NTuple{5,Any}, ::Val{2}) = - args[3] * args[5] * (args[2]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillAR), args::NTuple{5,Any}, ::Val{3}) = (args[1]^args[5])   /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])
Symbolics.derivative(::typeof(HillAR), args::NTuple{5,Any}, ::Val{4}) = - args[3] * args[5] * (args[3]^(args[5]-1)) * (args[1]^args[5])  /  (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
Symbolics.derivative(::typeof(HillAR), args::NTuple{5,Any}, ::Val{5}) = args[3] * (args[1]^args[5])  *  (log(args[1])*(args[2]^args[5] + args[4]^args[5]) - (args[2]^args[5])*log(args[2]) - (args[4]^args[5])*log(args[4]))   /   (args[1]^args[5] + args[2]^args[5] + args[4]^args[5])^2
