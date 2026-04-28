using EnumClasses, Printf

abstract type AbstractUpdate end
struct TMIN <:AbstractUpdate end
struct TMAX <:AbstractUpdate end
struct TMAX_SPAN <:AbstractUpdate end
struct TMIN_SPAN <:AbstractUpdate end

@enumclass Update tmin tmax tmax_span tmin_span

mutable struct Ranges{T<:AbstractUpdate}
    tmin
    tmax
    is
    ie
    step
    function Ranges(tmin,tmax,step,update::Update.Type)
        update == Update.tmin && return new{TMIN}(tmin,tmax,tmin,tmax,step)
        update == Update.tmax && return new{TMAX}(tmin,tmax,tmin,tmax,step)
        if update == Update.tmax_span
            length(tmax) >1 || error("Expected more than one tmax")
            return new{TMAX_SPAN}(tmin,tmax,tmin,tmax[1],step)
        end
        if update == Update.tmin_span
            length(tmin)>1 || error("Expected more than one tmin")
            return new{TMIN_SPAN}(tmin,tmax,tmin[1],tmax,step)
        end
    end
    Ranges(tmin,tmax,update::Update.Type) =  Ranges(tmin,tmax,1,update)
end

function Base.iterate(r::Ranges{T}) where T<:AbstractUpdate
    r.is=r.tmin[1]
    r.ie=r.tmax[1]
    return range(r.is, r.ie, step=r.step), 1
end

function Base.iterate(r::Ranges{TMIN},state)
    if r.ie-r.is == 1
        r.ie -=1
        r.ie == r.tmin && return nothing
        r.is = r.tmin
        return range(r.is,r.ie,step = r.step), state+1
    end
    r.is +=1
    return range(r.is,r.ie,step=r.step), state+1
end

function Base.iterate(r::Ranges{TMAX},state)
    if r.ie-r.is == 1
        r.is +=1
        r.is == r.tmax && return nothing
        r.ie = r.tmax
        return range(r.is,r.ie,step = r.step), state+1
    end
    r.ie -=1
    return range(r.is,r.ie,step=r.step), state+1
end

function Base.iterate(r::Ranges{TMAX_SPAN},state)
    if r.ie-r.is ==1
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
    if r.ie-r.is ==1
        i = findlast(r.is.== r.tmin)
        i == lastindex(r.tmin) && return nothing
        r.is = r.tmin[i+1]
        r.ie = r.tmax
        return range(r.is,r.ie,step=r.step),state+1
    end
    r.ie -=1
    return range(r.is,r.ie,step=r.step),state+1
end


Base.IteratorSize(::Ranges{T}) where T <:AbstractUpdate = Base.SizeUnknown()
Base.IteratorSize(::Ranges{TMIN})  = Base.HasLength()
Base.IteratorSize(::Ranges{TMAX})  = Base.HasLength()
Base.IteratorSize(::Ranges{TMAX_SPAN}) = Base.HasLength()

function Base.length(r::Ranges{T}) where T<: Union{TMIN,TMAX}
    L = r.tmax-r.tmin+1
    return div(L*(L-1),2)
end

Base.length(r::Ranges{TMAX_SPAN}) =  sum(t-r.tmin for t in r.tmax)
Base.length(r::Ranges{TMIN_SPAN}) =  sum(r.tmax-t for t in r.tmin)

function fit_scan(model,xdata,ydata,npar;
                  tmin = 1, tmax = lastindex(xdata),
                  update = Update.tmin,
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
    ranges = Ranges(tmin,tmax,update)
    fits =Vector{Union{FitRes,Missing}}(missing, length(ranges))
    for (idx,r) in enumerate(ranges)
        length(r)<=npar && continue;

        println(logfile, "="^80, "\n")
        Printf.@printf logfile "Fit range  [%.5e, %.5e]\n" xdata[r[1]] xdata[r[end]]
        fit = fit_routine(model,xdata[r],ydata[r],npar, W=W,C=C,
                          guess=guess,logfile=logfile, corr=corr)
        fits[idx] = FitRes(;fit...,range = (xdata[r[1]],xdata[r[end]]))
    end
    flush(logfile)
    return collect(skipmissing(fits))
end
