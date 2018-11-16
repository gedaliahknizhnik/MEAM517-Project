function [dw] =  chaplyginSleigh(t,w,params)
% The state w = [x,y,th,ph,dx,dy,dth,dph,Tapp]

% Select control input
TappVals = params.Tapp;

for ii = 1:size(TappVals,1)
    if (TappVals(ii,1) > t)
        break;
    end
end

dt = TappVals(ii,1)-TappVals(ii-1,1);
dTapp = TappVals(ii,2)-TappVals(ii-1,2);
Tapp = TappVals(ii-1,2) + (dTapp)*(t-TappVals(ii-1,1))/dt;

% Calculate dynamics
dw = f(w,Tapp,params);

waitbar(t/params.maxTime,params.waitBar,'Solving ODE')

end