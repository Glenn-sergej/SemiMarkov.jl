function set_semimarkov(Settings, StateTransProb, StateTrans, t, dst)

    C = zeros(Settings["trans"], length(t), length(t));
    dC = zeros(Settings["trans"], length(t), length(t));
    F = zeros(Settings["states"], length(t), length(t));
    for nt in 1:Settings["trans"]
        for st in 1:length(t)
            for rt in 1:length(t)
                if dst[nt,1] == "Exponential"
                    dC[nt,st,rt] = StateTransProb[rt,nt]*pdf(Exponential(dst[nt,2]),t[st])
                    C[nt,st,rt] = StateTransProb[rt,nt]*cdf(Exponential(dst[nt,2]),t[st])
                elseif dst[nt,1] == "Weibull"
                    dC[nt,st,rt] = StateTransProb[rt,nt]*pdf(Weibull(dst[nt,2], dst[nt,3]),t[st])
                    C[nt,st,rt] = StateTransProb[rt,nt]*cdf(Weibull(dst[nt,2], dst[nt,3]),t[st])
                elseif dst[nt,1] == "Lognormal"
                    dC[nt,st,rt] = StateTransProb[rt,nt]*pdf(Lognormal(dst[nt,2], dst[nt,3]),t[st])
                    C[nt,st,rt] = StateTransProb[rt,nt]*cdf(Lognormal(dst[nt,2], dst[nt,3]),t[st])
                end
            end
        end
        from = StateTrans[nt,1];
        F[from, :, :] = F[from, :, :] + C[nt, :, :]
    end
    return C, dC, F
end

function replace_nan(x)
    for i = eachindex(x)
        if isnan(x[i])
            x[i] = 0
        end
    end
end
        