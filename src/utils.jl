function TIC(f::FitRes{NT}) where NT <:Union{NamedTuple,Nothing}
    return f.chi2 - 2f.chiexp
end

function AIC(f::FitRes{NT},ndata::Int64) where NT<:Union{NamedTuple,Nothing}
    return f.chi2 + 2length(f.par) - 2ndata
end
