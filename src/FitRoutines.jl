module FitRoutines

using ADerrors, LsqFit, LinearAlgebra, ForwardDiff


include("types.jl")
include("pvalue.jl")
include("fit_functions.jl")
include("find_fit.jl")
include("fit_scan.jl")
include("utils.jl")

export pvalue, fit_routine, FitRes, fit_scan
end # module FitRoutines
