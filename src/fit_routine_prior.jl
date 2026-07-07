@doc raw"""

y = model(x,p,prior)

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
    return fit
end
