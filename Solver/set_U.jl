function set_U(Settings, t, dC, StateTrans)
    U = zeros(Settings["states"]*length(t), Settings["states"]*length(t))
    ID = Matrix(1.0I,Settings["states"], Settings["states"])
    for rt in 1:length(t)
        w = set_w(rt)
        for st in 1:rt            
            Ψ = set_Psi(Settings, st, rt, dC, StateTrans)
            i = (rt-1)*Settings["states"]+1:rt*Settings["states"]
            j = (st-1)*Settings["states"]+1:st*Settings["states"]
            if rt==st
                U[i,j] = ID-w[st]*Settings["dT"]*Ψ
            else
                U[i,j] = -w[st]*Settings["dT"]*Ψ
            end
        end
    end
    return U
end

function set_Psi(Settings, st, rt, dC, StateTrans)
    Ψ = zeros(Settings["states"], Settings["states"])
    for x in 1:Settings["trans"]
        from = StateTrans[x,1]
        to = StateTrans[x,2]
        Ψ[to, from] = dC[x,st,rt]
    end
    return Ψ
end

function set_w(x)
    w = zeros(x)
    for i in 1:length(w)
        if i % 2 == 0
           w[i] = 4/3 
        else
            w[i] = 2/3
    end end
    w[1] = 1/3
    w[end] = 1/3

    return w
end






