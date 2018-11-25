function [K,u] = findK(t, w, ts, Ls, z_states, z_controls, R, polys, params)

    [w_star,u_star] = getDes(t,ts,z_states,z_controls,polys);

    [A,B] = linSys(w_star,u_star,params);
 
    L = reshape(deval(Ls,t),params.nx,params.nx);
    S = L*L';
    K = R^(-1)*B'*S;
    
    u = u_star + K*(w_star-w);
end