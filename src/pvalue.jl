@doc """
    __Q(nu,chi,nmc)

Internal Function, it computes the pvalue given the matrix `nu`, the `chi2` and the number of random sample `nmc`

"""
function __Q(nu,chi,nmc)
    Q = 0.0
    N,_ = size(nu)
    z = zeros(N)
    z2 = similar(z)
    eig = abs.(LinearAlgebra.eigvals(nu))
    eps = 1e-14 * maximum(eig)
    for i in eachindex(eig)
        eig[i]  = eig[i]>eps ? eig[i] : 0.0
    end
    @inbounds for i in 1:nmc
        for j in eachindex(z)
            z[j] = randn()
            z2[j] = z[j]^2
        end
        aux = sum(z2[i]*eig[i] for i in eachindex(eig))
        if aux < chi
            Q+=1.0
        end
    end
    return 1.0 - Q/nmc
end

@doc raw"""
    pvalue(chisq::Function,
    chi2::Float64,
    xp::AbstractVector{Float64},
    data::AbstractVector{uwreal},
    W::Union{Vector{Float64},Array{Float64,2}} = Vector{Float64}();
    wpm = Dict{Int64,Vector{Float64}}()
    nmc::Int64 = 5000,
    C::AbstractMatrix{Float64} = Matrix{Float64}(undef,0,0))

Computes the p-value of a previously done fit, using as input the `\chi^2` observed from the fit, the fit parameters and the fitted data.
The p-value for a given `\chi^2` is the probability of, given the data you have, finding such a `\chi^2` or worse from a fit, and still
have the data well described by the fit function. `nmc` is the number of MC samples used to estimate the p-value integral, default is 5000.
By now it only works with a vector for weights (containing the diagonal of W)

```@example
function fit_defs(f::Function,x,W)
    chisq(p,d)=sum((d-f(x,p)).^2 .*W)
    return chisq
end

@.fun(x,p) = p[1] + p[2] * x
chisq = fit_defs(fun, x, 1.0 ./ err.(y) .^ 2)
fit = curve_fit(fun, x, value.(y), 1.0 ./ err.(y) .^ 2, [0.5, 0.5])
(up, chiexp) = fit_error(chisq, coef(fit), y)

wpm = Dict{Int64,Vector{Float64}}()
wpm[1] = [-1.0,-1.0,4-0,-1.0]
Q = pvalue(chisq, chi2, value.(up), y, wpm; W = 1.0 ./ err.(y) .^ 2, nmc=10000)
#Q = pvalue(chisq, chi2, value.(up), y; W = 1.0 ./ err.(y) .^ 2, nmc=10000)
#Q = pvalue(chisq, chi2, value.(up), y)
```
"""
function pvalue(chisq::Function,
    chi2::Float64,
    xp::AbstractVector{Float64},
    data::AbstractVector{uwreal},
    W::AbstractVector{Float64}=Vector{Float64}();
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} =  Dict{Int64,Vector{Float64}}(),
    nmc::Int64 = 5000,
    C::AbstractMatrix{Float64} = Matrix{Float64}(undef,0,0))

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data
    if (m-n<0)
      error("more parameters than data")
    end

    xav = zeros(Float64, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end
    ccsq(x::Vector) = chisq(view(x, 1:n), view(x, n+1:n+m))
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{4}());
    end

    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)

    if length(W) ==0
      W = zeros(Float64, m)
      for i in 1:m
        if (data[i].err == 0.0)
          uwerr(data[i], wpm)
          if (data[i].err == 0.0)
            error("Zero error in fit data")
          end
        end
        W[i] = 1.0 / data[i].err^2
      end
    end

    n = size(hess, 1) - m
    hm = view(hess, 1:n, n+1:n+m)
    sm = Array{Float64, 2}(undef, n, m)

    for i in 1:n, j in 1:m
        sm[i,j] = hm[i,j] / sqrt.(W[j])
    end

    maux = Array{Float64,2}(undef,n,n)

    LinearAlgebra.mul!(maux,sm, sm')

    hi = LinearAlgebra.pinv(maux)
    Px = let
        _Px = zeros(Float64,m,m)
        LinearAlgebra.mul!(sm,hi,hm)
        LinearAlgebra.mul!(_Px,-hm',sm)
        for i in 1:m
            _Px[i,i] = W[i] + _Px[i,i]
        end
        _Px
    end

    length(C)!=0 || (C=ADerrors.cov_sym(data))

    nu = let
        aux = sqrt(C)
        aux2 = similar(aux)
        _nu = similar(aux)
        LinearAlgebra.mul!(aux2,aux,Px)
        LinearAlgebra.mul!(_nu,aux2,aux)
        _nu
    end

    return __Q(nu,chi2,nmc)
end


function pvalue(chisq::Function,
    chi2::Float64,
    xp::AbstractVector{Float64},
    data::AbstractVector{uwreal},
    W::AbstractArray{Float64,2};
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} =  Dict{Int64,Vector{Float64}}(),
    nmc::Int64 = 5000,
    C::AbstractMatrix{Float64} = GammaMethod.cov_AD(data,wpm))

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data
    if (m<n)
        error("More parameters that data in pvalue")
    end

    xav = zeros(Float64, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end

     ccsq(x::Vector) = chisq(view(x, 1:n), view(x, n+1:n+m))
    if (n+m < 4)
        cfg =  ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{1}());
    else
        cfg =  ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{4}());
    end


    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)

    m = length(data)
    n = size(hess, 1) - m

    Lm = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(W))
    Li = LinearAlgebra.inv(Lm.L)

    hm = view(hess, 1:n, n+1:n+m)
    sm = Array{Float64,2}(undef,n,m)
    LinearAlgebra.mul!(sm,hm,Li')

    maux = Array{Float64,2}(undef,n,n)
    LinearAlgebra.mul!(maux,sm,sm')

    hi   = LinearAlgebra.pinv(maux)
    Px = let
        _Px = Array{Float64,2}(undef,m,m)
        LinearAlgebra.mul!(sm,hi,hm)
        LinearAlgebra.mul!(_Px,hm',sm)
        for i in eachindex(_Px)
            _Px[i] = W[i] - _Px[i]
        end
        _Px
    end

    nu = let
        aux=sqrt(C)
        aux2 = similar(aux)
        _nu = similar(aux)
        LinearAlgebra.mul!(aux2,aux,Px)
        LinearAlgebra.mul!(_nu,aux2,aux)
        _nu
    end
    return __Q(nu,chi2,nmc)
end
