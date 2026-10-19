__get_scanned_region(xdata, idx_min, idx_max) = xdata[minimum(idx_min)], xdata[maximum(idx_max)]
function __get_scanned_region(xdata, idx_min::Function, idx_max)
    idx = maximum(idx_max)
    xdata[idx_min(idx)], xdata[idx]
end
function __get_scanned_region(xdata, idx_min, idx_max::Function)
    idx = maximum(idx_min)
    xdata[idx], xdata[idx_max(idx)]
end

@doc raw"""
     find_first(model::Function, xdata::AbstractArray, ydata::AbstractArray, npar::Int64; condition::Function, idx_min=1, idx_max=length(xdata), step = 1, W = Floatt64[], logfile = NoPrint(), wpm = Dict{Int64,Vector{Float64}}(), C = Matrix{Float64}(undef,0,0), guess = fill(0.5,npar), corr=false)

Find the first fit which fulfills `condition`. It fits the data in `ydata` with the function `model` by calling the method [`fit_routine`](@ref). The fits are performed over the ranges `Ranges(idx_min,idx_max,step)`, see also [`Ranges`](@ref).
The function `condition` accepts a [`FitRes`](@ref) object and return `true` if the fit is successful.
The fit are performed in the order specify by `Ranges(idx_min,idx_max,step)` and the first successful fit is returned. An optimal choice of `idx_min`, `idx_max` and `step`, might find the desired fit faster.

A log of the fits can be produced by passing the file handle in `logfile`. By default, no logfile is printed.
If an error occurs during the exectution of [`fit_routine`](@ref), the  error is catched, a error message is printed in terminal and in the logfile, then the process continue with the next fit.

## Return type
if a fit is found, it return a  `FitRes` object with `idx_range` and `x_range` storing the range  of the successful fit. `idx_range` stores the range in index notation, `x_range` stores the range in `xdata`, such that `x_range = xdata[idx_range]`.
If no fit is found, it return `nothing`
"""
function find_first(model::Function, xdata::AbstractArray, ydata::AbstractArray, npar::Int64;
                    condition::Function,
                    idx_min = 1, idx_max = lastindex(xdata),
                    step=1,
                    logfile = FitRoutines.NoPrint(),
                    wpm   = Dict{Int64,Vector{Float64}}(),
                    guess = fill(0.5, npar),
                    corr = false,
                    Ctol = -1.0,
                    Ctol_method = :default,
                    cov_func = ADerrors.cov_sym)

    println(logfile,"Starting fits: find first mode")
    let (xmin, xmax) = __get_scanned_region(xdata,idx_min,idx_max)
        Printf.@printf logfile "Scanning region xdata in [%.12e, %.12e]\n" xmin xmax
    end
    ranges = Ranges(idx_min,idx_max,step)
    println(logfile, ranges)
    flush(logfile)
    for r in ranges
        length(r)<=npar && continue;
        println(logfile, "="^80, "\n")
        Printf.@printf logfile "Fit range [%.5e, %.5e]\n" xdata[r[1]] xdata[r[end]]
        fit = try
            fit_routine(model, xdata[r], ydata[r], npar, wpm=wpm, guess=guess, logfile=logfile, corr=corr, Ctol = Ctol, Ctol_method = Ctol_method, cov_func = cov_func)
        catch e
            e isa Base.InterruptException && rethrow(e)
            Printf.@printf logfile "Error during fit routine: \n"
            showerror(logfile, e,catch_backtrace())
            println(logfile)
            @debug Printf.@printf "Error in fit range [%.5e,%.5e]\n" xdata[r[1]] xdata[r[end]]
            showerror(stdout, e)
            println()
            continue;
        end
        condition(fit) && return FitRes(;fit..., x_range = xdata[r[[1,end]]], idx_range = r[[1,end]])
        flush(logfile)
    end
    return nothing
end

@doc raw"""
     find_first(model::Function, xdata::AbstractArray, ydata::AbstractArray, npar::Int64; condition::Function, idx_min=1, idx_max=length(xdata), step = 1, W = Floatt64[], logfile = NoPrint(), wpm = Dict{Int64,Vector{Float64}}(), C = Matrix{Float64}(undef,0,0), guess = fill(0.5,npar), corr=false)

Find the first fit which fulfills `condition`. It fits the data in `ydata` with the function `model` by calling the method [`fit_routine`](@ref). The fits are performed over the ranges `Ranges(idx_min,idx_max,step)`, see also [`Ranges`](@ref).
The function `condition` accepts a [`FitRes`](@ref) object and return `true` if the fit is successful.
The fit are performed in the order specify by `Ranges(idx_min,idx_max,step)` and the first successful fit is returned. An optimal choice of `idx_min`, `idx_max` and `step`, might find the desired fit faster.

A log of the fits can be produced by passing the file handle in `logfile`. By default, no logfile is printed.
If an error occurs during the exectution of [`fit_routine`](@ref), the  error is catched, a error message is printed in terminal and in the logfile, then the process continue with the next fit.

## Return type
If a fit is found, it return a `Vector{FitRes}`, each  with `idx_range` and `x_range` storing the range  of the successful fit. `idx_range` stores the range in index notation, `x_range` stores the range in `xdata`, such that `x_range = xdata[idx_range]`.
If no fit is found, it return `nothing`
"""
function find_all(model::Function,
                  xdata::AbstractArray,
                  ydata::AbstractArray,
                  npar::Int64;
                  condition::Function,
                  idx_min = 1,
                  idx_max = lastindex(xdata),
                  logfile = FitRoutines.NoPrint(),
                  wpm = Dict{Int64,Vector{Float64}}(),
                  guess = fill(0.5, npar),
                  Ctol = -1.0,
                  Ctol_method = :default,
                  cov_func = ADerrors.cov_sym)
    println(logfile,"Starting fits: find first mode")
    let (xmin,xmax)  = __get_scanned_region(xdata,idx_min,idx_max)
        Printf.@printf logfile "Scanning region xdata in [%.12e, %.12e]\n" xmin xmax
    end
    ranges = Ranges(idx_min,idx_max)
    println(logfile, ranges)
    flush(logfile)
    res = FitRes[];
    for r in ranges
        length(r)<=npar && continue;
        println(logfile, "="^80, "\n")
        Printf.@printf logfile "Fit range  [%.5e, %.5e]\n" xdata[r[1]] xdata[r[end]]
        fit = try
            fit_routine(model, xdata[r], ydata[r], npar, wpm=wpm, guess=guess, logfile=logfile, corr=corr, Ctol = Ctol,
                        Ctol_method = Ctol_method, cov_func = cov_func)
        catch e
            Printf.@printf logfile "Error during fit routine: %s\n" e.msg
            @debug Printf.@printf "Error in fit range [%.5e,%.5e]\n" xdata[r[1]] xdata[r[end]]
            continue;
        end
        condition(fit) && push!(res,FitRes(;fit..., x_range = xdata[r[[1,end]]], idx_range = r[[1,end]]))
        flush(logfile)
    end
    return length(res) == 0 ? nothing : res
end

export find_first, find_all
