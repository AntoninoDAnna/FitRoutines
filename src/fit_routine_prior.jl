@doc raw"""
     fit_routine_prior(model::Function,xdata)
"""
function fit_routine_prior(model::Function,xdata::AbstractArray{<:Real},ydata::AbstractArray{uwreal},prior::AbstractArray{uwreal},npar::Int64;
                           wpm::Union{Dict{Int64,Vector{Float64}}, Dict{String,Vector{Float64}}} = Dict{Int64,Vector{Float64}}(),
                           corr::Bool = true,
                           W::AbstractArray{Float64}=Vector{Float64}(),
                           guess::AbstractVector{Float64} = fill(0.5,npar),
                           logfile = NoPrint(),
                           C::AbstractArray{Float64} = Matrix{Float64}(undef,0,0))
    nobs = length(ydata)
    nprior = length(prior)
    length(xdata)==nobs || error("xdata and ydata musth be of the same size")
    function F(i::Vector,par)
        res = Vector{Any}(undef,nobs+nprior)
        res[1:nobs] .= model(xdata,par[1:npar],par[npar+1:npar+nprior])
        res[nobs+1:nobs+nprior] .= par[npar+1:npar+nprior]
        return res
    end
    guess = vcat(guess,value.(prior))
    y = vcat(ydata,prior)
    fit = fit_routine(F, collect(eachindex(y)), y, length(guess),
                      guess=guess, W=W, C=C, corr=corr,
                      logfile=NoPrint(),
                      wpm=wpm)
    return FitRes(;fit...,prior=prior)
end

function fit_routine_prior(model::Function,
                           xdata::AbstractArray{<:Real},
                           ydata::AbstractArray{uwreal},
                           prior::AbstractArray{<:Tuple{<:Real,<:Real}},
                           npar::Int64;
                           wpm::Union{Dict{Int64,Vector{Float64}}, Dict{String,Vector{Float64}}} = Dict{Int64,Vector{Float64}}(),
                           corr::Bool = true,
                           W::AbstractArray{Float64}=Vector{Float64}(),
                           guess::AbstractVector{Float64} = fill(0.5,npar),
                           logfile = NoPrint(),
                           C::AbstractArray{Float64} = Matrix{Float64}(undef,0,0))
    nobs = length(ydata)
    nprior = length(prior)
    length(xdata)==nobs || error("xdata and ydata musth be of the same size")
    function F(i::Vector,par)
        res = Vector{Any}(undef,nobs+nprior)
        res[1:nobs] .= model(xdata,par[1:npar],par[npar+1:npar+nprior])
        res[nobs+1:nobs+nprior] .= par[npar+1:npar+nprior]
        return res
    end
    vprior = getfield.(prior,1)
    eprior = getfield.(prior,2)
    guess = vcat(guess,vprior)
    y = vcat(value.(ydata),vprior)
    [uwerr(y, wpm) for y in ydata]
    if length(W) ==0
        if corr
            W = zeros(ndata+nprior,ndata+nprior)
            C = length(C) ==0 ? ADerrors.cov_sym(ydata) : C
            W[1:ndata,1:ndata] = LinearAlgebra.pinv(C);
            for i in ndata+1:ndata+prior
                W[i,i] = 1/eprior[i]^2
            end
            W = (W'+W)*0.5
        else
            W = vcat([1/y.err^2 for y in ydata],1 ./eprior.^2)
        end
    elseif size(W,1) == ndata
        if corr
            _W = zeros(ndata+nprior,ndata+nprior)
            _W[1:ndata,1:ndata] .= W
            for i in ndata+1:ndata+prior
                _W[i,i] = 1/eprior[i]^2
            end
            W = _W
        else
            W = vcat(W,1 ./eprior^2)
        end
    end
    fit = LsqFit.curve_fit(F,xdata,y,W,guess)

    chisq(p,d) = get_chi(F,xdata,p,vcat(d,vprior),W)
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

    return FitRes(par=par,chi2=chi2,chiexp=chiexp,pval=pval,prior=prior)
end

function fit_routine_prior(models::AbstractArray{<:Function},
                           xdata::AbstractArray{<:AbstractArray},
                           ydata::AbstractArray{<:AbstractArray{ADerrors.uwreal}},
                           prior::AbstractArray,
                           npar::Int64;
                           W = Float64[],
                           C = Matrix{Float64}(undef,0,0),
                           wpm::AbstractDict = Dict{Int64,Vector{Float64}}(),
                           corr::Bool = true,
                           guess::AbstractArray = fill(0.5,npar),
                           logfile = NoPrint(),
                           )
    Nmodel = length(models)
    ndata = length.(xdata)
    Ntotal = sum(ndata);
    if Nmodel != length(xdata) != length(ydata)
        error("[FitRoutines] Error: you need the same number of model, xdata and ydata")
    end
    if !(length(W) == 0 || length(W)==Nmodel || length(W) == Ntotal^2)
        error("[FitRoutines] Error: You need to pass a W matrix for model, for all the data or none")
    end
    if !(length(C) == 0 || length(C)==Nmodel || length(C) == Ntotal^2)
        error("[FitRoutines] Error: You need to pass a C matrix for model, for all the data or none")
    end
    if !all(length.(ydata).== ndata)
        error("[FitRoutines] Error: Mismatch in xdata and ydata. Make sure that the data corresponds")
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

    function Model(x::Vector,p::AbstractVector,prior::AbstractVector;models =models, Nmodel::Int64=Nmodel,ndata::Vector{Int64}=ndata)
        res = Vector{Any}(undef,Nmodel)
        res[1] = models[1](x[1:ndata[1]],p,prior)
        is= 1
        ie= 0
        for i in 1:Nmodel
            ie += ndata[i]
            rd = is:ie
            res[i] = models[i](x[rd],p,prior)
            is = ie +1
        end
        return vcat(res...)
    end

    fit = fit_routine_prior(Model,X,Y,prior,npar,guess = guess,W=WW,C=CC,corr=corr,logfile=NoPrint(),wpm=wpm)
    compute_error_for_logging(logfile, fit.par)
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
