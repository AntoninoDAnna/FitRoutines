@doc raw"""
    fit_scan(model,xdata,ydata,npar;
             tmin=1, tmax=lastindex(xdata), kwargs...)

It exectutes all the fit compatible with `tmin`, `tmax`. Return a  `Vector{FitRes}
"""
function fit_scan(model,xdata,ydata,npar;
                  tmin = 1, tmax = lastindex(xdata),
                  W = Float64[],
                  logfile = NoPrint(),
                  wpm   = Dict{Int64,Vector{Float64}}(),
                  C     = Matrix{Float64}(undef, 0,0),
                  guess = fill(0.5, npar),
                  corr = false)

    println(logfile,"Starting Fit Scan")
    Printf.@printf logfile "Scanning region [%.12e, %.12e]\n" xdata[tmin[1]] xdata[tmax[1]]
    ranges = Ranges(tmin,tmax)
    println(logfile, ranges)
    flush(logfile)
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
        fits[idx] = FitRes(;fit...,x_range = [xdata[r[1]],xdata[r[end]]], idx_range = r[[1,end]])
    end
    flush(logfile)
    return collect(skipmissing(fits))
end
