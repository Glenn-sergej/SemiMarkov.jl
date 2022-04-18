include("set_u.jl")

function set_p(Settings, rt, h, F, Init, ns)

    p = 0;
    #Initial value
    p = p + Init[ns]*(1-F[ns, rt, 1]);
    w = set_w(rt);
    for st in 1:rt
        p = p + Settings["dT"]*w[st]*h[st,ns]*(1-F[ns, rt, st]);
    end
    return p
end