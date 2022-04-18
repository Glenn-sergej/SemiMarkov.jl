function set_U(Settings, t, dC, StateTrans)
    U = zeros(Settings["states"]*length(t), Settings["states"]*length(t))
    ID = Matrix(I,Settings["states"], Settings["states"])
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

function set_w(t)
    weights = zeros(t)
    if t==1
        weights = [0];
    elseif t==2
        weights = [1/2, 1/2];
    elseif t==3
        weights = [1/3, 4/3, 1/3];
    elseif t==4
        weights = [3/8, 9/8, 9/8, 3/8];
    elseif t==5
        weights = (2/45)*[7, 32, 12, 32, 7];
    elseif t==6
        weights = (5/288)*[19,75,50,50,75,19];
    elseif t==7 
        weights = (1/140)*[41,216,27,272,27,216,41];
    elseif t==8
        weights = (7/17280)*[751,3577,1323,2989,2989,1323,3577,751];
    else
        weights = ones(t)*48;
        weights[1:4,1] = [17;59;43;49];
        weights[end-3:end,1] = [49;43;59;17];
        weights = (1/48).*weights;
    end
    return weights
end







