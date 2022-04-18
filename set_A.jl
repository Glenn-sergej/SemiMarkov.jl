function set_A(Settings, StateTrans, Init, dC, t)
    A = zeros(Settings["states"]*length(t));
    for x in 1:length(t)
        for tr in 1:Settings["trans"]
            from = StateTrans[tr,1];
            to = StateTrans[tr,2];
            adj_pos = Settings["states"]*(x-1) + to;
            A[adj_pos] = A[adj_pos] + Init[from]*dC[tr,x,1];
        end
    end
    return A
end