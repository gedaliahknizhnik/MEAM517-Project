x_0 = [0,0,0,0,0,0,0,0];
x_f = [30,0,0,0,0,0,0,0];
x_f = [-35.1616, -1.1016, 0, 0, 0, 0, 0, 0]

ind1 = [1,1,1,1,1,1,1,1];
ind2 = [1,1,0,0,0,0,0,0];

% Setup simulation parameters.
params.m     = 10;
params.B     = 1e1;
params.C     = 1e1;
params.a     = 1;
params.Tresb = 0;
params.Tresc = 0;
params.Fx    = 0;
params.Fy    = 0;
params.Fxw   = 0;
params.Fyw   = 0;
params.maxTime = 10;