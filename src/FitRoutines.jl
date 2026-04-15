module FitRoutines

using ADerrors, LsqFit, LinearAlgebra


include("types.jl")
include("pvalue.jl")
include("fit_functions.jl")

export pvalue, fit_routine
end # module FitRoutines
