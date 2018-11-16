%% Run the simulation of the Chaplygin Sleigh
% Still figuring out the best way to pass the forces to the ODE function.

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

% Set up ODE solver

initth  = 0;
initVel = 0;
initdx  = initVel*cos(initth);
initdy  = initVel*sin(initth);
inits   = [0,0,initth,0,initdx,initdy,0,0];
tInt    = [0,params.maxTime];


%% Solve the ODE

params.Tapp = [0,   -7.0195
   1, -5.8956
   2, -4.2119
   3, -2.5271
   4, -0.8424
   5,  0.8424
   6,  2.5271
   7,  4.2119
   8,  5.8956
   9,  7.0195];

W = waitbar(0,'Solving ODE'); params.waitBar = W;
[t,z] = ode45(@chaplyginSleigh,tInt,inits,[],params);
close(W)

%% Validate the solution by checking the constraint:
thT = z(:,3); dxT = z(:,5); dyT = z(:,6); dthT = z(:,7); 
constraint = dxT.*sin(thT) - dyT.*cos(thT) + params.a*dthT;

disp('Max Constraint Error:')
max(constraint)

%% Plot the outputs

close all
figure, hold on
plot(z(:,1),z(:,2),'Linewidth',2)
title('Trajectory','Interpreter','Latex')
xlabel('$x$ (m)','Interpreter','Latex')
ylabel('$y$ (m)','Interpreter','Latex')
axis equal

figure, hold on
plot(t, z(:,1));
plot(t, z(:,2));
plot(t, z(:,3));
plot(t, z(:,4));
title('Positions','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value (m,rad)','Interpreter','Latex')
legend({'$x$','$y$','$\theta$','$\phi$'},'Interpreter','Latex')

figure, hold on
plot(t, z(:,5));
plot(t, z(:,6));
plot(t, z(:,7));
plot(t, z(:,8));
title('Velocities','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value (m/s,rad/s)','Interpreter','Latex')
legend({'$\dot{x}$','$\dot{y}$','$\dot{\theta}$','$\dot{\phi}$'},'Interpreter','Latex')
