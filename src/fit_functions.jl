@doc """
    get_chi(model,x,par,data,W::AbstractMatrix)
    get_chi(model,x,par,data,W::AbstractVector)

It computes the chi square of a fit. `model` is assumed to be a function called as `f(x[i],par)`
If `W::AbstractVector` the chi square is computed for an uncorrelated fit
If `W::AbstractMatrix` the chi square is computed for a correlated fit

It assumes that `x`, `data` and `W` have compatible dimensions.
"""
function get_chi(model,x,par,data,W::AbstractMatrix)
    return sum((data[i]-model(x[i],par)) * W[i,j] * (data[j]-model(x[j],par)) for i in axes(W,1), j in axes(W,2))
end

function get_chi(model,x,par,data,W::AbstractVector)
    return sum((data[i]-model(x[i],par))^2 * W[i] for i in eachindex(W))
end

@doc raw"""
        fit_routine(model::Function,xdata::AbstractArray{<:Real},ydata::AbstractArray{ADerrors.uwreal},npar; kwargs...)
        fit_routine(model::Function,xdata::AbstractArray{uwreal},ydata::AbstractArray{ADerrors.uwreal},npar; kwargs...)

It executes a fit of the `ydata` with the fit function `model`.

# Arguments:

  - `model`: Fit function used to fit the data. the function assumes that its called as `model(xdata,parameters)`.
  - `xdata`: indipendent variable in the fit. It can be an `AbstractArray{<:Real}`, that is the "default" method, or
             a `AbstractArray{uwreal}`. In the latter case, the method builds the function
``math
F(q_i) = \begin{cases}
         model(q_{npar+i},q[1:npar])  \quad \text{if} i<= length(xdata)\\
         q_i \quad \text{elsewhere}
       \end{cases}
``
             and then call the default method with dummy `xdata` and `vcat(ydata,xdata)` as new `ydata`
  - `ydata`: Data variable that carries the error
  - `npar`: number of parameters used in the fit.

# Keyword parameters:

  - `wpm`: window parameter used by ADerrors to computes the errors. Default to `Dict{Int64,Vector{Float64}}()`.
           This is used both when computing the error of the data and when computing the covariance matrix through `ADerrors.cov`
           It accept a dictionary with ensemble id as keys (either `String` or `Int64`) and `Vector{Float64}` as value.
           See also ADerrors documentation.

  - `W::AbstractVecOrMat{Float64}`: Weight matrix/vector. If `W` is a `Vector` the fit is perfomed uncorrelated ,otherwise is correlated.
                                    default to `Vector{Float64}()`

  - `corr::Bool`: flag used to generated `W` when not given. Default to `true`. If `false`, `W = [1/y.err^2 for y in ydata]`.
                  If `true` the covariance matrix of the data is used to compute `W` as `W = C^{-1}`.

  - `C::AbstractMatrix`: covariance matrix used to compute W. Default to `Matrix{Float64}(undef, 0, 0)`.
                         If `W` is not given, and `corr` is `true`, then it is used to compute `W`.
                         If need to computed `W` but not given,  then `C = ADerrors.cov(ydata,wpm)`. See also ADerrors documentation.
                         If available, it will be used to better estimate the fit pvalue and the expected chi-square,
                         but it will not be computed if not available at this point.

  - `guess::Vector{Float64}`: Initial guess for the parameter used in the fit routine. If not given, default to `fill(0.5,npar)`.
                              If `xdata isa Vector{uwreal}`, the mean value of `xdata` are used as guess for the relevant `q_i` parameters.

  - `logfile`: handle to a logfile. If `nothing`, no log will be written. It can be `IO` file or a custom log structure that has
               an overloading for `print()` and `println()`.

# Returns
It returns a `NamedTuple` with names:
  - `:par`:    fit parameter as `Vector{uwreal}`
  - `:chi2`:   chisquare,
  - `:chiexp`: chisquare expected
  - `:pval`:   pvalue
"""
function fit_routine(model::Function,
                     xdata::AbstractArray{<:Real},
                     ydata::AbstractArray{uwreal},
                     npar::Int64;
                     wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}}=Dict{Int64,Vector{Float64}}(),
                     corr::Bool = true,
                     W::AbstractArray{Float64}=Vector{Float64}(),
                     guess::AbstractVector{Float64} = fill(0.5,npar),
                     logfile = NoPrint(),
                     C::AbstractArray{Float64} = Matrix{Float64}(undef,0,0));

    nobs = length(ydata)
    if length(xdata)!=nobs
        error("xdata and ydata must be the of same size")
    end

    [uwerr(y, wpm) for y in ydata]
    if length(W) ==0
        if corr
            C = length(C) ==0 ? ADerrors.cov_sym(ydata) : C
            W = LinearAlgebra.pinv(C);
            W = (W'+W)*0.5
        else
            W = [1/y.err^2 for y in ydata]
        end
    end

    fit = LsqFit.curve_fit(model,xdata,ADerrors.value.(ydata),W,guess)

    chisq(p,d) = get_chi(model,xdata,p,d,W)
    chi2 = sum(fit.resid.^2)
    par = fit_error(chisq,LsqFit.coef(fit),ydata,wpm,W,false)
    # if the fit is correlated  we have already C, no need to recompute it.
    chiexp  = ADerrors.chiexp(chisq,
                              LsqFit.coef(fit),
                              ydata,
                              wpm,
                              C=length(C)==0 ? nothing : C,
                              W=W)

    pval = pvalue(chisq,chi2,LsqFit.coef(fit),ydata,W,wpm=wpm,C=C)

    compute_error_for_logging(logfile,par);
    println(logfile, "fit_routine output:")
    println(logfile, "\t\tFit routine in [$(xdata[1]),...,$(xdata[end])]")
    println(logfile, "\t\tparameters :")
    for p in par
        println(logfile,"\t\t\t",p)
    end
    println(logfile, "\t\tchi2:       ", chi2)
    println(logfile, "\t\tchiexp:     ", chiexp)
    println(logfile, "\t\tpvalue:     ", pval)

    return FitRes(par=par,chi2=chi2,chiexp=chiexp,pval=pval)
end

"""
        check_size(M,n)

Check if `M` has sizes all equal to n or if is empty
"""
function check_sizes(M,n)
    if length(M) ==0 || all(size(M) .== n)
        return true
    end
    return false
end

function fit_routine(model::Function,
                     xdata::AbstractArray{ADerrors.uwreal},
                     ydata::AbstractArray{ADerrors.uwreal},
                     npar::Int64;
                     W::AbstractArray = Float64[],
                     wpm::AbstractDict = Dict{Int64,Vector{Float64}}(),
                     corr::Bool = true,
                     guess::AbstractArray = fill(0.5,npar),
                     logfile = NoPrint(),
                     C::AbstractArray = Float64[])
    Nx = length(xdata)
    if (Nx != length(ydata))
        error("[juobs] Error: xdata and ydata must have the same dimension")
    end
    Ndata = 2Nx; ## in total we have 2*lengt(xdata) data points
    data = vcat(ydata,xdata)
    if !check_sizes(W,Ndata)
        error("[juobs] Error: W does not have the correct sizes")
    end
    if !check_sizes(C,Ndata)
        error("[juobs] Error: C does not have the correct sizes")
    end

    if corr
        C = length(C) ==0 ? GammaMethod.cov_AD(data) : C
        W = length(W) ==0 ? LinearAlgebra.pinv(W) : W
    else
        W = length(W) ==0 ? [1/yy.err^2 for yy in data] : W
    end

    f(q,p) = [i<=Nx ? model(q[i],p) : q[i-Nx] for i in 1:Ndata]
    fit_func(x,p) = f(view(p,npar+1:npar+Nx),view(p,1:npar))
    new_guess = [guess; value.(xdata)]
    fit = fit_routine(fit_func,zeros(Ndata),data,npar+Nx,
                      guess=new_guess, wpm=wpm,corr=corr,
                      W=W,C=C,logfile=nothing)

    compute_error_for_logging(logfile,par);
    println(logfile, "fit_routine output:")
    println(logfile, "\t\txdata carries uncertainty")
    println(logfile, "\t\tFit routine in [$(xdata[1].mean),...,$(xdata[end].mean)]")
    println(logfile, "\t\tFit parameters:")
    for i in 1:npar
        println(logfile,"\t\t\t",fit.par[i])
    end
    println(logfile, "\t\tauxiliary fit parameter:")
    println(logfile, "\t\t\txdata[i] \t q[i]")
    for i in 1:Nx
        println(logfile, "\t\t\t",xdata[i]," \t ", fit.par[npar+i])
    end
    println(logfile, "\t\tchi2:       ", fit.chi2)
    println(logfile, "\t\tchiexp:     ", fit.chiexp)
    println(logfile, "\t\tpvalue:     ", fit.pval)
    return fit
end

@doc raw"""
        fit_routine(model::AbstractArray{Function},
                    xdata::AbstractArray{<:AbstractArray},
                    ydata::AbstractArray{<:AbstractArray{ADerrors.uwreal}},npar; kwargs...)


It executes a combined fit of the `ydata` with the fit functions in `model`. The method define the general model `F(X_{mi},p) = model[m](xdata[i],p)`
and call `fit_routine(F,vcat(xdata),vcat(ydata),npar;kwargs....)`. The `kwargs` parameters are accordingly updated to reflect the new function.

# Returns
It returns a `NamedTuple` with names:
  - `:par`:    fit parameter as `Vector{uwreal}`
  - `:chi2`:   chisquare,
  - `:chiexp`: chisquare expected
  - `:pval`:   pvalue
"""
function fit_routine(models::AbstractArray{Function},
                     xdata::AbstractArray{<:AbstractArray},
                     ydata::AbstractArray{<:AbstractArray{ADerrors.uwreal}},
                     npar::Int64;
                     W = Float64[],
                     C = Float64[],
                     wpm::AbstractDict = Dict{Int64,Vector{Float64}}(),
                     corr::Bool = true,
                     guess::AbstractArray = fill(0.5,npar),
                     logfile = NoPrint(),
                     )
    Nmodel = length(models)
    ndata = length.(xdata)
    Ntotal = sum(ndata);
    if Nmodel != length(xdata) != length(ydata)
        error("[juobs] Error: you need the same number of model, xdata and ydata")
    end
    if !(length(W) == 0 || length(W)==Nmodel || length(W) == Ntotal^2)
        error("[juobs] Error: You need to pass a W matrix for model, for all the data or none")
    end
    if !(length(C) == 0 || length(C)==Nmodel || length(C) == Ntotal^2)
        error("[juobs] Error: You need to pass a C matrix for model, for all the data or none")
    end
    if !all(length.(ydata).== ndata)
        error("[juobs] Error: Mismatch in xdata and ydata. Make sure that the data corresponds")
    end
    X = vcat(xdata...)
    Y = vcat(ydata...)

    function make_matrix(M::Vector{<:AbstractMatrix})
        aux = zeros(Ntotal,Ntotal)
        ie = 0
        for m in 1:Nmodel
            is = ie+1
            ie += ndata[m];
            aux[is:ie,is:ie] .= M[m]
        end
        return aux
    end
    make_matrix(M::T where {T<:AbstractMatrix}) = M

    WW = length(W)==0 ? W : make_matrix(W);
    CC = length(C)==0 ? C : make_matrix(C);

    function Model(x,p)
        res = Vector{Any}(undef,Nmodel)
        res[1] = models[1](x[1:ndata[1]],p)
        is= 1
        ie= 0
        for i in 1:Nmodel
            ie += ndata[i]
            rd = is:ie
            res[i] = models[i](x[rd],p)
            is = ie +1
        end
        return vcat(res...)
    end

    fit = fit_routine(Model,X,Y,npar,guess = guess,W=WW,C=CC,corr=corr,logfile=nothing,wpm=wpm)


    compute_error_for_logging(logfile,par)
    flag = typeof(xdata) <: AbstractArray{<:AbstractArray{ADerrors.uwreal}}
    println(logfile, "fit_routine output:")
    if flag
        println(logfile, "\t\txdata carries uncertainty")
    end
    println(logfile, "\t\tGlobal fit with $Nmodel functions")
    for nm in 1:Nmodel
        println(logfile,"---Model $(nm)---")
        println(logfile, "\t\tFit routine in [$(xdata[nm][1]),...,$(xdata[nm][end])]")
    end
    println(logfile, "\t\tFit parameters:")
    for i in 1:npar
        println(logfile,"\t\t\t",fit.par[i])
    end
    if flag
        off = npar
        println(logfile, "\t\tauxiliary fit parameter:")
        for nm in 1:Nmodel
            println(logfile,"---Model $(nm)---")
            println(logfile, "\t\t\txdata[m][i] \t q[i]")
            for i in eachindex(xdata[nm])
                println(logfile, "\t\t\t",xdata[nm][i]," \t ", fit.par[off+i])
            end
            off += npar+length(xdata[nm])
        end
    end
    println(logfile, "\t\tchi2:       ", fit.chi2)
    println(logfile, "\t\tchiexp:     ", fit.chiexp)
    println(logfile, "\t\tpvalue:     ", fit.pval)
    flush(logfile)
    return fit
end
