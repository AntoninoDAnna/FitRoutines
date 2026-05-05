using EnumClasses, Printf

abstract type AbstractUpdate end
struct ALL <:AbstractUpdate end
struct TMAX_SPAN <:AbstractUpdate end
struct TMIN_SPAN <:AbstractUpdate end
struct TMINMAX_SPAN <:AbstractUpdate end

mutable struct Ranges{T<:AbstractUpdate}
    tmin
    tmax
    is
    ie
    step
    Ranges(tmin, tmax, step=1) = Ranges(tmin,tmax,step,eltype(tmin),eltype(tmax))
    Ranges(tmin::T,tmax::U,step,::Type{T},::Type{U}) where {T<:Number,U<:Number} = new{ALL}(tmin,tmax,tmin,tmax,step)
    Ranges(tmin::T,tmax::U,step,::Type{S},::Type{V}) where {T,U,S<:Number,V<:Number} = new{TMINMAX_SPAN}(tmin,tmax,tmin[1],tmax[1],step)
    Ranges(tmin::T,tmax::U,step,::Type{T},::Type{V}) where {T<:Number,U,V<:Number}  = new{TMAX_SPAN}(tmin,tmax,tmin,tmax[1],step)
    Ranges(tmin::T,tmax::U,step,::Type{V},::Type{U}) where {T, U<:Number,V<:Number} = new{TMIN_SPAN}(tmin,tmax,tmin[1],tmax,step)
end

function Base.iterate(r::Ranges{T}) where T<:AbstractUpdate
    r.is=r.tmin[1]
    r.ie=r.tmax[1]
    return range(r.is, r.ie, step=r.step), 1
end

function Base.iterate(r::Ranges{ALL},state)
    if r.ie-r.is <= 1
        r.ie -=1
        r.ie == r.tmin && return nothing
        r.is = r.tmin
        return range(r.is,r.ie,step = r.step), state+1
    end
    r.is +=1
    return range(r.is,r.ie,step=r.step), state+1
end

function Base.iterate(r::Ranges{TMAX_SPAN},state)
    if r.ie-r.is <=1
        i = findlast(r.ie.== r.tmax)
        i == lastindex(r.tmax) && return nothing
        r.is = r.tmin
        r.ie = r.tmax[i+1]
        return range(r.is,r.ie,step=r.step),state+1
    end
    r.is +=1
    return range(r.is,r.ie,step=r.step),state+1
end

function Base.iterate(r::Ranges{TMIN_SPAN},state)
    if r.ie-r.is <=1
        i = findlast(r.is.== r.tmin)
        i == lastindex(r.tmin) && return nothing
        r.is = r.tmin[i+1]
        r.ie = r.tmax
        return range(r.is,r.ie,step=r.step),state+1
    end
    r.ie -=1
    return range(r.is,r.ie,step=r.step),state+1
end

function Base.iterate(r::Ranges{TMINMAX_SPAN},state)
    if r.is == r.tmin[end]
        r.ie == r.tmax[end] && return nothing
        i= findlast(r.ie==t for t in r.tmax)
        r.ie=r.tmax[i+1]
        r.is=r.tmin[1]
    else
        i= findlast(r.is==t for t in r.tmin)
        r.is = r.tmin[i+1]
    end
    return range(r.is,r.ie,step=r.step), state
end

Base.IteratorSize(::Ranges{T}) where T <:AbstractUpdate = Base.SizeUnknown()


@doc raw"""
    fit_scan(model,xdata,ydata,npar;
             tmin=1, tmax=lastindex(xdata), kwargs...)

It exectutes all the fit compatible with `tmin`, `tmax`. Return a  `Vector{FitRes}
"""
function fit_scan(model,xdata,ydata,npar;
                  tmin = 1, tmax = lastindex(xdata),
                  W = Float64[],
                  logfile = FitRoutines.NoPrint(),
                  wpm   = Dict{Int64,Vector{Float64}}(),
                  C     = Matrix{Float64}(undef, 0,0),
                  guess = fill(0.5, npar),
                  corr = false)

    println(logfile,"Starting Fit Scan")
    Printf.@printf logfile "Scanning region [%.12e, %.12e]\n" xdata[1] xdata[end]
    Printf.@printf logfile "Scanning update %s\n" string(update)

    vdata = Vector{NamedTuple}()  ## vector of fit results
    ranges = Ranges(tmin,tmax)
    fits =Vector{Union{FitRes,Missing}}(missing, length(ranges))
    for (idx,r) in enumerate(ranges)
        length(r)<=npar && continue;

        println(logfile, "="^80, "\n")
        Printf.@printf logfile "Fit range  [%.5e, %.5e]\n" xdata[r[1]] xdata[r[end]]
        fit = try
            fit_routine(model,xdata[r],ydata[r],npar, W=W,C=C,
                        guess=guess,logfile=logfile, corr=corr)
        catch e
            Printf.@printf logfile "Error during fit routine: %s\n" e.msg
            @error "Error in fit range [%.5e,%.5e]\n" xdata[r[1]] xdata[r[end]]
            continue;
        end
        fits[idx] = FitRes(;fit...,range = (xdata[r[1]],xdata[r[end]]))
    end
    flush(logfile)
    return collect(skipmissing(fits))
end
