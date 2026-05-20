function find_first(model,xdata,ydata,npar;
                    condition,
                    tmin = 1, tmax = lastindex(xdata),
                    W = Float64[],
                    logfile = FitRoutines.NoPrint(),
                    wpm   = Dict{Int64,Vector{Float64}}(),
                    C     = Matrix{Float64}(undef, 0,0),
                    guess = fill(0.5, npar),
                    corr = false)
    println(logfile,"Starting fits: find first mode")
    Printf.@printf logfile "Scanning region [%.12e, %.12e]\n" xdata[tmin[1]] xdata[tmax[1]]
    ranges = Ranges(tmin,tmax)
    println(logfile, ranges)
    flush(logfile)
    for r in ranges
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
        condition(fit) && return FitRes(;fit..., range = (xdata[r[1]], xdata[r[end]]))
        flush(logfile)
    end
    return nothing
end

function find_all(model,xdata,ydata,npar;
                    condition,
                    tmin = 1, tmax = lastindex(xdata),
                    W = Float64[],
                    logfile = FitRoutines.NoPrint(),
                    wpm   = Dict{Int64,Vector{Float64}}(),
                    C     = Matrix{Float64}(undef, 0,0),
                    guess = fill(0.5, npar),
                    corr = false)
    println(logfile,"Starting fits: find first mode")
    Printf.@printf logfile "Scanning region [%.12e, %.12e]\n" xdata[tmin[1]] xdata[tmax[1]]
    ranges = Ranges(tmin,tmax)
    println(logfile, ranges)
    flush(logfile)
    res = FitRes[];
    for r in reverse(ranges)
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
        condition(fit) && push!(res,FitRes(;fit..., range = (xdata[r[1]], xdata[r[end]])))
        flush(logfile)
    end
    return res
end
