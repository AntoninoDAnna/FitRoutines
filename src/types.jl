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
