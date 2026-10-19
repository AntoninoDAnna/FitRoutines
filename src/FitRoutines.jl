module FitRoutines

using ADerrors, LsqFit, LinearAlgebra, ForwardDiff


include("types.jl")
include("Gamma_method.jl")
include("fit_functions.jl")
include("fit_routine_prior.jl")
include("find_fit.jl")
include("fit_scan.jl")
include("utils.jl")

export fit_routine, FitRes, fit_scan
end # module FitRoutines
