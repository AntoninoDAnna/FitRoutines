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

import Base.getproperty

function Base.getproperty(v::FitRes{<:NamedTuple{names,types}}, name::Symbol) where {names, types}
    (name in names) && return getfield(v.extra, name)
    getfield(v,name)
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

import Printf:format
format(::NoPrint,f,x...) = nothing

## auxiliary function for uwreal

compute_error_for_logging(logfile, p) = uwerr.(p)
compute_error_for_logging(logfile::NoPrint,p) = nothing


## Ranges
using Printf

abstract type AbstractUpdate end
struct ALL <:AbstractUpdate end
struct TMAX_SPAN <:AbstractUpdate end
struct TMIN_SPAN <:AbstractUpdate end
struct TMINMAX_SPAN <:AbstractUpdate end

mutable struct Ranges{T<:AbstractUpdate}
    tmin
    tmax
    is
    ie
    step
    Ranges(tmin, tmax, step=1) = Ranges(tmin,tmax,step,eltype(tmin),eltype(tmax))
    Ranges(tmin::T,tmax::U,step,::Type{T},::Type{U}) where {T<:Number,U<:Number} = new{ALL}(tmin,tmax,tmin,tmax,step)
    Ranges(tmin::T,tmax::U,step,::Type{S},::Type{V}) where {T,U,S<:Number,V<:Number} = new{TMINMAX_SPAN}(tmin,tmax,tmin[1],tmax[1],step)
    Ranges(tmin::T,tmax::U,step,::Type{T},::Type{V}) where {T<:Number,U,V<:Number}  = new{TMAX_SPAN}(tmin,tmax,tmin,tmax[1],step)
    Ranges(tmin::T,tmax::U,step,::Type{V},::Type{U}) where {T, U<:Number,V<:Number} = new{TMIN_SPAN}(tmin,tmax,tmin[1],tmax,step)
end

function Base.iterate(r::Ranges{T}) where T<:AbstractUpdate
    r.is=r.tmin[1]
    r.ie=r.tmax[1]
    return range(r.is, r.ie, step=r.step), 1
end

function Base.iterate(r::Ranges{ALL},state)
    if r.ie-r.is <= 1
        r.ie -=1
        r.ie == r.tmin && return nothing
        r.is = r.tmin
        return range(r.is,r.ie,step = r.step), state+1
    end
    r.is +=1
    return range(r.is,r.ie,step=r.step), state+1
end

function Base.iterate(r::Ranges{TMAX_SPAN},state)
    if r.ie-r.is <=1
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
    if r.ie-r.is <=1
        i = findlast(r.is.== r.tmin)
        i == lastindex(r.tmin) && return nothing
        r.is = r.tmin[i+1]
        r.ie = r.tmax
        return range(r.is,r.ie,step=r.step),state+1
    end
    r.ie -=1
    return range(r.is,r.ie,step=r.step),state+1
end

function Base.iterate(r::Ranges{TMINMAX_SPAN},state)
    if r.is == r.tmin[end]
        r.ie == r.tmax[end] && return nothing
        i= findlast(r.ie==t for t in r.tmax)
        r.ie=r.tmax[i+1]
        r.is=r.tmin[1]
    else
        i= findlast(r.is==t for t in r.tmin)
        r.is = r.tmin[i+1]
    end
    return range(r.is,r.ie,step=r.step), state
end

Base.IteratorSize(::Ranges{T}) where T <:AbstractUpdate = Base.HasSize()

function Base.length(r::Ranges{ALL})
    L = r.tmax-r.tmin+1
    div(L*(L-1),2)
end

Base.length(r::Ranges{TMIN_SPAN}) = sum(r.tmax<=t ? 0 : r.tmax-t for t in r.tmin)
Base.length(r::Ranges{TMAX_SPAN}) = sum(t<=r.tmin ? 0 : t-r.tmin for t in r.tmax)
Base.length(r::Ranges{TMINMAX_SPAN}) = sum(T<=t ? 0 : T-t for T in r.tmax, t in r.tmin )

Base.print(io::IO,r::Ranges{TMIN_SPAN}) = @printf io "Tmin span in [%d,%d], tmax start in %d" r.tmin[1] r.tmin[end] r.tmax

Base.print(io::IO,r::Ranges{TMAX_SPAN}) = @printf io "Tmax span in [%d,%d], tmin starts in %d" r.tmax[1] r.tmax[end] r.tmin

Base.print(io::IO,r::Ranges{TMINMAX_SPAN}) = @printf io "Tmin span in [%d,%d], tmax span in [%d,%d]" r.tmin[1] r.tmin[end] r.tmax[1] r.tmax[end]

Base.print(io::IO,r::Ranges{ALL}) = @printf io "All fit in [%d,%d]" r.tmin r.tmax
