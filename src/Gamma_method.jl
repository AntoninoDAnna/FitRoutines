import ADerrors, ForwardDiff, LinearAlgebra

@doc raw"""
      calc_gamma(a1::ADerrors.uwreal [,a2::Aderrors.uwreal] ,id::Union{Int64,String}[, ws::ADerrors.wspace])

Computes the Gamma function associated to the ensemble `id` without the bias included by ADerrors for t>=0

If `id` is not associated to a Montecarlo chain, it return a Vector with the correlation between `a1` and `a2`
If at least one obs does not depend on `id`, it return `[0.0]`
"""

function calc_gamma(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::Int64,ws::ADerrors.wspace)
    idx  = ws.map_ids[id]
    nd   = ws.fluc[idx].nd
    if !(id in a1.ids) || !(id in a2.ids)
        return [0.0]
    end

    if nd==1
        function get_v(a,id=id,ws=ws)
            v = 0
            for i in eachindex(a.prop)
                if (a.prop[i] && (ws.map_nob[i] == id))
                    v+=a.der[i].*ws.fluc[i].delta[1]
                end
            end
            return v
        end
        return [get_v(a1)*get_v(a2)]
    end

    nd_eff = div(nd, ws.fluc[idx].ibn)
    (nt,_) = findmax(ws.fluc[idx].ivrep)
    nt = div(nt, 2*ws.fluc[idx].ibn)
    nrep = length(ws.fluc[idx].ivrep)
    ftemp1 = Dict{Int64,Vector{Complex{Float64}}}()
    ftemp2 = Dict{Int64,Vector{Complex{Float64}}}()

    for i in 1:nrep
        ns = length(ws.fluc[idx].fourier[i])
        ftemp1[i] = zeros(Complex{Float64}, ns)
        ftemp2[i] = zeros(Complex{Float64}, ns)
    end

    gamma = zeros(Float64, nt)

    for i in eachindex(a1.prop)
        if (a1.prop[i] && (ws.map_nob[i] == id))
            for k in 1:nrep
                ftemp1[k] = ftemp1[k] + a1.der[i]*ws.fluc[i].fourier[k]
            end
        end
    end

    for i in eachindex(a2.prop)
        if (a2.prop[i] && (ws.map_nob[i] == id))
            for k in 1:nrep
                ftemp2[k] = ftemp2[k] + a2.der[i]*ws.fluc[i].fourier[k]
            end
        end
    end

    for k in 1:nrep
        ftemp2[k] .= ftemp1[k].*conj(ftemp2[k])
        ADerrors.FFTW.ifft!(ftemp2[k])

        for ig in 1:min(nt,length(ftemp2[k]))
            gamma[ig] = gamma[ig] + real(ftemp2[k][ig])
        end
    end

    for ig in 1:nt
        nrcnt = count(map(x -> div(x, 2*ws.fluc[idx].ibn), ws.fluc[idx].ivrep) .> ig-1)
        gamma[ig] = gamma[ig] / (nd_eff - nrcnt*(ig-1))
    end
    return gamma
end

calc_gamma(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::String,ws::ADerrors.wspace) = calc_gamma(a1,a2,ws.str2id[id],ws)
calc_gamma(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::String) = calc_gamma(a1,a2,ADerrors.wsg.str2id[id],ADerrors.wsg)
calc_gamma(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::Int64)  = calc_gamma(a1,a2,id,ADerrors.wsg)
calc_gamma(a::ADerrors.uwreal,id::Int64,ws::ADerrors.wspace)   = calc_gamma(a,a,id,ws)
calc_gamma(a::ADerrors.uwreal,id::String,ws::ADerrors.wspace)  = calc_gamma(a,a,ws.str2id[id],ws)
calc_gamma(a::ADerrors.uwreal,id::String) = calc_gamma(a,a,ADerrors.wsg.str2id[id],ADerrors.wsg)
calc_gamma(a::ADerrors.uwreal,id::Int64)  = calc_gamma(a,a,id,ADerrors.wsg)

function drho(rho::Vector{Float64},nd::Int64,iw::Int64=0)
    iw = iw ==0 ? ADerrors.wopt_ulli(nd,ADerrors.DEFAULT_STAU,rho) : iw
    nt = length(rho);
    dr = zeros(nt);
    for t in 1:nt
        is = max(1, t-iw-2) + 1
        ie = is + iw-1
        for m in is:ie
            aux=0.0
            if m <= nt
                aux = -2*rho[t]*rho[m]
            end
            if abs(m-t)<nt
                aux += rho[abs(m-t)+1]
            end
            if m+t-1<=nt
                aux += rho[m+t-1]
            end
            dr[t] +=aux^2
        end
        dr[t] = sqrt(dr[t]/nd)
    end
    return dr
end

"""
      bias(a1::ADerrors.uwreal [,a2::ADerrors.uwreal], id::Union{Int64,String} [, ws::Aderrors.wspace]; iw::Int64=0)

  it return the bias associated to `id` that ADErrors includes in the estimate of the autocorrelation function

      B = Γ(0) + 2*∑_{t=0}^{iw}Γ(t)

  to retrieve the unbias error of `a`

      σ_a^{unbias}  = √(a.err^2 - (2W-1)/N^2 B)

  if `iw = 0` the integration window is determined by the automatic windowing procedure with S_τ = 4.0

  """
function bias(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::Int64,ws::ADerrors.wspace;iw::Int64=0)
    idx  = ws.map_ids[id]
    nd   = ws.fluc[idx].nd
    if nd ==1
        return 0.0
    end
    gamma = calc_gamma(a1,a2,id,ws)
    return gamma[1] + 2.0*sum(gamma[2:iw])
end

bias(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::String,ws::ADerrors.wspace; iw::Int64=0)= bias(a1,a2,ws.str2id[id],ws,iw=iw)
bias(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::Int64; iw::Int64=0) = bias(a1,a2,id,wsg,iw=iw)
bias(a1::ADerrors.uwreal,a2::ADerrors.uwreal,id::String; iw::Int64=0)= bias(a1,a2,ADerrors.wsg.str2id[id],wsg,iw=iw)

function bias(a::ADerrors.uwreal,id::Int64,ws::ADerrors.wspace;iw::Int64=0)
    idx  = ws.map_ids[id]
    nd   = ws.fluc[idx].nd
    if nd ==1
        return 0.0
    end
    gamma = calc_gamma(a,id,ws)
    iw = iw ==0 ?  ADerrors.wopt_ulli(nd,ADerrors.DEFAULT_STAU,gamma) : iw
    return gamma[1] + 2.0*sum(gamma[2:iw])
end

bias(a::ADerrors.uwreal,id::String,ws::ADerrors.wspace; iw::Int64=0) = bias(a,ws.str2id[id],ws,iw=iw)
bias(a::ADerrors.uwreal,id::Int64; iw::Int64=0)  = bias(a,id,ADerrors.wsg,iw=iw)
bias(a::ADerrors.uwreal,id::String; iw::Int64=0) = bias(a,ADerrors.wsg.str2id[id],ADerrors.wsg,iw=iw)

function bias(vobs::Vector{ADerrors.uwreal}, id::Int64, ws::ADerrors.wspace;)
    b = zeros(length(vobs),length(vobs))
    iw = maximum([ADerrors.window(obs,id) for obs in vobs])
    for i in 1:length(vobs)-1, j in i+1:length(vobs)
        b[i,j] = bias(vobs[i],vobs[j],id,ws,iw=iw)
        b[j,i] = b[i,j]
    end
    for i in eachindex(vobs)
        b[i,i] = bias(vobs[i],id,ws,iw = 0)
    end
    return b
end

bias(vobs::Vector{ADerrors.uwreal}, id::String, ws::ADerrors.wspace) = bias(vobs,ws.str2id[id],ws)
bias(vobs::Vector{ADerrors.uwreal}, id::Int64)  = bias(vobs,id,ADerrors.wsg)
bias(vobs::Vector{ADerrors.uwreal}, id::String) = bias(vobs,ADerrors.wsg.str2id[id],ADerrors.wsg)


@doc raw"""
      cov_pe(vobs::Vectir{ADerrors.uwreal}[,ws:ADerrors.wspace])

  Computes the unbiased covariance matrix of `vobs` à la `pyerrors`.

  It first computes the naive correlation matrix and then it multiplies it with the errors
  of the `vobs`. This method is exact if

      tau_int^ij = sqrt(tau_int^i tau_int^j)

  It assumes that `ADerrors.uwerr` has been run on `vobs`
  This method is stable in the sense that the covarince of two `ADerrors.uwreal` depends only
  on them and not on the others `ADerrors.uwreal` that are in `vobs`

  See also [`cov_min`](@ref) [`cov_ad`](@ref)
  """
function cov_pe(vobs::Vector{ADerrors.uwreal}, ws::ADerrors.wspace = ADerrors.wsg)
    ids = ADerrors.unique_ids_multi(vobs, ws)
    nid = length(ids)
    nobs = length(vobs)
    cov = zeros(Float64, nobs,nobs)

    # naive covarince matrix
    for j in 1:nid
        idx  = ws.map_ids[ids[j]]
        nd_eff = div(ws.fluc[idx].nd, ws.fluc[idx].ibn)
        for n1 in 1:nobs, n2 in n1:nobs
            gamma = calc_gamma(vobs[n1],vobs[n2],ids[j],ws)
            cov[n1,n2] +=gamma[1]/nd_eff
        end
    end
    for n1 in 1:nobs-1, n2 in n1+1:nobs ## symmetrizing
        cov[n2,n1] = cov[n1,n2]
    end
    corr = [cov[i,j]/sqrt(cov[i,i]*cov[j,j])  for i in 1:nobs, j in 1:nobs]
    cov  = [corr[i,j]*vobs[i].err*vobs[j].err for i in 1:nobs, j in 1:nobs]
    return cov
end

cov_pe(vobs::Vector{ADerrors.uwreal}, wpm::AbstractDict, ws::ADerrors.wspace=ADerrors.wsg) = cov_pe(vobs, ws)


@doc raw"""
      cov_min(vobs::Vector{ADerrors.uwrela} [ws::ADerrors.wspace])

  It computes the unbiased covariance matrix of `vobs` using a mixed method.

  For each id, it first determines the smallest integration window used in the data
  and computes "à la Aderrors" the correlation matrix, then it computes the covarince
  matrix "à la pyerrors, that is, it multipy the covarince matrix with the errors of `vobs`

  It assumes that `ADerrors.uwerr` has been run on `vobs`
  This method is not stable, in the sense that the covarince between two `ADerrors.uwreal`
  depend also on the other `ADerrors.uwreal` included in `vobs`

  See also [`cov_pe`](@ref) [`cov_ad`](@ref)
  """
function cov_min(vobs::Vector{ADerrors.uwreal},ws::ADerrors.wspace)
    ids = ADerrors.unique_ids_multi(vobs, ws)
    nid = length(ids)
    iw  = fill(typemax(Int64),nid)
    for obs in vobs, j in eachindex(obs.ids), i in eachindex(ids)
        if (obs.ids[j] == ids[i])
            iw[i] = min(iw[i], obs.cfd[j].iw)
        end
    end

    nobs = length(vobs)
    cov = zeros(Float64, nobs,nobs)

    for j in 1:nid
        idx  = ws.map_ids[ids[j]]
        nd_eff = div(ws.fluc[idx].nd, ws.fluc[idx].ibn)

        for n in eachindex(vobs)
            gamma = calc_gamma(vobs[n],ids[j],ws)
            cov[n,n] += length(gamma)==1 ? gamma[1]/nd_eff : (gamma[1] + 2.0*sum(gamma[2:iw[j]]))/nd_eff
        end
        for n1 in 1:nobs-1, n2 in n1+1:nobs
            gamma12 = calc_gamma(vobs[n1],vobs[n2],ids[j],ws)
            if length(gamma12) == 1
                cov[n1,n2] +=gamma12[1]/nd_eff
                continue;
            end
            gamma21 = calc_gamma(vobs[n2],vobs[n1],ids[j],ws)
            cov[n1,n2] += (gamma12[1] + sum(gamma12[2:iw[j]]) + sum(gamma21[2:iw[j]]))/nd_eff
        end
    end

    for n1 in 1:nobs-1, n2 in n1+1:nobs
        cov[n2,n1] = cov[n1,n2]
    end

    corr = [cov[i,j]/sqrt(cov[i,i]*cov[j,j]) for i in 1:nobs, j in 1:nobs]
    cov = [corr[i,j]*vobs[i].err*vobs[j].err for i in 1:nobs, j in 1:nobs]
    return cov
end

cov_min(vobs::Vector{ADerrors.uwreal}) = cov_min(vobs,ADerrors.wsg)
