struct FitRes{NT<:Union{NamedTuple,Nothing}}
    par::Union{Vector{uwreal},Missing}
    pval::Float64
    chi2::Float64
    chiexp::Float64
    extra::NT
    function FitRes(;par=[], pval = 0.0, chi2 = 0.0, chiexp = 0.0, extra...)
        isempty(extra) && return new{Nothing}(par,pval,chi2,chiexp,nothing)
        E = (;extra...)
        new{typeof(E)}(par,pval,chi2,chiexp,(;extra...))
    end
end

import Base.iterate

function Base.iterate(f::FitRes{NT}, state=1) where NT<:NamedTuple
    state == 1 && (return (:par => f.par), state+1)
    state == 2 && (return (:pval => f.pval),state+1)
    state == 3 && (return (:chi2 => f.chi2),state+1)
    state == 4 && (return (:chiexp => f.chiexp),state+1)
    idx = state-4
    idx > length(A) && return nothing
    return (keys(fit.extra)[i] => fit.extra[i]),state+1
end

function Base.iterate(f::FitRes{Nothing}, state=1)
    state == 1 && (return (:par => f.par), state+1)
    state == 2 && (return (:pval => f.pval),state+1)
    state == 3 && (return (:chi2 => f.chi2),state+1)
    state == 4 && (return (:chiexp => f.chiexp),state+1)
    return nothing
end


import Base: print
function __print_minimum_fit__(io::IO,f::FitRes{NT}) where NT <:Union{Nothing,NamedTuple}
    if ismissing(f.par)
        print(io, "Failed fit: No parameter is available")
        return
    end
    println(io, "Fit with $length(f.par) parameters: ")
    for (i,p) in enumerate(f.par)
        println(io,"\tp_",i-1,"=\t",p)
    end
    println(io,"\tpval =\t",f.pval)
    println(io,"\tchi^2=\t",f.chi2,"\tchi^2_exp=",f.chiexp)
    print(io,"\tchi^2/chi_exp  = ",f.chi2/f.chiexp)
end

print(io::IO, f::FitRes{Nothing}) = __print_minimum_fit__(io,f)

function print(io::IO,f::FitRes{NT}) where NT<:NamedTuple
    __print_minimum_fit__(io,f)
    ismissing(f.par) && return
    for k in keys(f.extra)
        println()
        print(io,"\t",k,"=\t",getfield(f.extra,k))
    end
    nothing
end
"""
        struct NoPrint

empty structure used to disable anything tied to a logging process.
"""
struct NoPrint end

import Base.print, Base.println, Base.flush
print(::NoPrint,x...) = nothing
println(::NoPrint,x...) = nothing
Base.flush(::NoPrint) = nothing

## auxiliary function for uwreal

compute_error_for_logging(logfile, p) = uwerr.(p)
compute_error_for_logging(logfile::NoPrint,p) = nothing
