function Soln = trajectoryLqr(t,xd,ud,linSys,Q,R,F,tol,params)
% Soln = trajectoryLqr(t,linSys,Q,R,F,tol)
%
% This function is used to solve the finite-horizon continuout-time linear
% quadratic regulator problem.
%
% INPUTS:
%   t = monotonically increasing vector of times at which you would like
%       the gain matrix to be calculated
%   linSys = function handle for time-varying linear dynamics
%       [A, B] = linSys(t)
%   Q = state cost matrix 
%   R = input cost matrix 
%   F = final state cost matrix 
%   tol = accuracy of riccati equation propagation
%
% OUTPUTS:
%   Soln = struct array with solution at each point in t
%   Soln(i).t = t(i);
%   Soln(i).K = gain matrix at t(i)
%   Soln(i).S = Cost to go
%   Soln(i).E = close-loop eigen-values for system at t(i)
%
% NOTES:
%
%   J = x'Fx + Integral {x'Qx + u'Ru} dt
%
% See Also LQR
%
% CREDIT to Matthew Kelly, modified by Gedaliah Knizhnik
% (https://www.mathworks.com/matlabcentral/fileexchange/54432-continuous-time-finite-horizon-lqr)

nState = size(Q,1);
nInput = size(R,1);

td = t;

userFun = @(t,z)rhs(t,td,xd,ud,z,linSys,Q,R,nState,params);
z0 = reshape(F,nState*nState,1);
tSpan = [t(end),t(1)];

options = odeset();
options.RelTol = tol;
options.AbsTol = tol;
sol = ode45(userFun,tSpan,z0);
z = deval(sol,t);

nSoln = length(t);
Soln(nSoln).t = 0;
Soln(nSoln).K = zeros(nState,nInput);
Soln(nSoln).S = zeros(nState,nState);
Soln(nSoln).E = zeros(nState,1);

for idx=1:nSoln
    i = nSoln-idx+1;
    zNow = z(:,i);
    tNow = t(i);
        xdNow = interp1(td,xd,tNow);
        udNow = interp1(td,ud,tNow);

    S = reshape(zNow,nState,nState);
    [A,B] = linSys(xdNow,udNow,params);
    K = R\(B'*S);
    Soln(i).t = tNow;
    Soln(i).K = K;
    Soln(i).S = S;
    Soln(i).E = eig(A-B*K);
end

end

function dz = rhs(t,td,xd,ud,z,linSys,Q,R,nState, params)
    xdNow = interp1(td,xd,t);
    udNow = interp1(td,ud,t);
    
    P = reshape(z,nState,nState);
    [A,B] = linSys(xdNow,udNow,params);
    dP = ricatti(A,B,Q,R,P);
    dz = reshape(dP,nState*nState,1);
end

function [dP, K] = ricatti(A,B,Q,R,P)
K = R\B'*P;
dP = -(A'*P + P*A - P*B*K + Q);
end


