module FitRoutines

using ADerrors, LsqFit, LinearAlgebra, ForwardDiff


include("types.jl")
include("pvalue.jl")
include("fit_functions.jl")
include("fit_scan.jl")


export pvalue, fit_routine, FitRes, fit_scan
end # module FitRoutines
