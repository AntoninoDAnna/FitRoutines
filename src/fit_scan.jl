using EnumClasses, Printf
@enumclass Update tmin tmax

mutable struct Ranges
    tmin::Int64
    tmax::Int64
    is::Int64
    ie::Int64
    step::Int64
    update::Update.Type
    Ranges(tmin,tmax,update) =  new(tmin,tmax,tmin,tmax,1,update)
    Ranges(tmin,tmax,step,update) =  new(tmin,tmax,tmin,tmax,step,update)
end

function Base.iterate(r::Ranges)
    r.is=r.tmin
    r.ie=r.tmax
    return range(r.is, r.ie, step=r.step), 1
end

function Base.iterate(r::Ranges,state)
    if r.ie-r.is == 1
        if r.update == Update.tmin
            r.ie -=1
            r.ie == r.tmin && return nothing
            r.is = r.tmin
        elseif r.update==Update.tmax
            r.is +=1
            r.is == r.tmax && return nothing
            r.ie = r.tmax
        end
        return range(r.is,r.ie,step = r.step), state+1
    end
    if r.update == Update.tmin
        r.is +=1
    elseif r.update==Update.tmax
        r.ie -=1
    end
    return range(r.is,r.ie,step=r.step), state+1
end

function Base.IteratorSize(r::Ranges)
    if r.update in (Update.tmin, Update.tmax)
        return Base.HasLength()
    end
    return Base.SizeUnknown()
end

function Base.length(r::Ranges)
    if r.update in (Update.tmin, Update.tmax)
        L = r.tmax-r.tmin+1
        return div(L*(L-1),2)
    end
    return Base.SizeUnknown()
end

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
