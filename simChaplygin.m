%% SIMULATE CHAPLYGIN BEANIE
% This function sets up a Chaplygin Beanie problem and solves it using
% trajectory optimization to set up a direct collocation problem and solve
% it using SNOPT. It then simulates the output using the exact equations of
% motion and plots the comparison.

close all

nx = 8;
nu = 1;
N = 31;
T = 30;
dt = T/(N-1);

w_0 = [0,0,0,0,0,0,0,0];
w_f = [5,0,0,0,0,0,0,0];

ind1 = [1,1,1,1,1,1,1,1];
ind2 = [1,1,1,0,0,0,1,1];

M = Inf;

% Setup simulation parameters.
params.nx      = nx;
params.nu      = nu;
params.m       = 2;
params.B       = 1e0;
params.C       = 1e0;
params.a       = 0.05;
params.Tresb   = 0;
params.Tresc   = 0;
params.Fx      = 0;
params.Fy      = 0;
params.Fxw     = 0;
params.Fyw     = 0;
params.maxTime = T;
params.k1      = 0.0;

%% Solve the Optimal Control Problem
% Use SNOPT to calculate the optimal trajectory.

tic
[z,F,INFO] = optimizeChaplygin(w_0,w_f,ind1,ind2,M,nx,nu,N,dt,params);
params.tSNOPT = toc;

%% Process Output
% Prepare data for graphing and further use.
num = size(z,1)/(nx+nu);

z_SNOPT = zeros(num,nx);
u_SNOPT = zeros(num,nu);

nx_indices = unwrap(1:nx);
nu_indices = nx + unwrap(1:nu);

for ii = 1:num
    x_indices = (ii-1)*(nx+nu) + nx_indices;
    u_indices = (ii-1)*(nx+nu) + nu_indices;
    z_SNOPT(ii,:)   = z(x_indices)';
    u_SNOPT(ii,:) = z(u_indices)';
end

t_SNOPT = zeros(size(z_SNOPT,1),1);
for ii=1:size(t_SNOPT)
    t_SNOPT(ii) = (ii-1)*dt;
end

% Find polynomial trajectories output by SNOPT (currently unused because
% gives bad outputs for LQR - maybe b/c there's an error in there).
polys = getPolys(z_SNOPT,u_SNOPT,nx,dt,params);

% max(F)

%% Simulate ODE with these inputs 
% No control applied at this point - differences in trajectory between this
% and the optimal trajectory are assumed to be a combination of:
% # Numerical integration errors
% # Dynamic infeasibilities allowed by SNOPT but rejected by the ODE.

% Calculate the torques in a way the ODE will accept (for linear
% interpolation).
Tapp = zeros(size(u_SNOPT,1),2);
for ii = 1:size(Tapp,1)
    Tapp(ii,:) = [(ii-1)*dt, u_SNOPT(ii)];
end

% Store the torques for passing to the ODE solver.
params.Tapp = Tapp;

% Solve the ODE
% W = waitbar(0,'Solving ODE'); params.waitBar = W;
    tic
[t_OL,z_OL] = ode45(@chaplyginSleigh,[0,T],w_0,[],params);
    params.tODE = toc;


%% Find an LQR controller for this optimal trajectory
% Idea is to minimize trajectory deviation due to dynamic infeasibilities,
% even if it means deviating in terms of the control.

% LQR Costs
Q = 10*eye(params.nx);
R = 0.1*eye(params.nu);

% Final Cost-to-Go (must be non-zero b/c current method inverts matrix).
F = 0.01*eye(size(Q));
F = reshape(F,params.nx^2,1);

% Solve Riccati equation by integrating backwards in time 
%   (Note that reshaping occurs inside the function since ODE45 solves 
%   vectors, not matrices).
Lfun = @(t,L) Ldot(t, L, Q, R, t_SNOPT, z_SNOPT, u_SNOPT, polys, params);
% W = waitbar(0,'Solving Riccati equation ODE'); params.waitBar = W;
    tic
Lsol = ode45(Lfun, [params.maxTime,0], F);
    params.tRiccati = toc;
%close(W)

% Solve the dynamics with the LQR controller included (should see better
% performance than simply using the controller from SNOPT).
LQRfun = @(t,w) dynamics(t,w, t_SNOPT, Lsol, z_SNOPT, u_SNOPT, R, polys, params);
    tic
[t_LQR,z_LQR] = ode45(LQRfun,[0,T],w_0);
    params.tLQR = toc;

% Find the controls at each step 
u_LQR = zeros(size(t_LQR));

for ii = 1:size(t_LQR)
    [~,u] = findK(t_LQR(ii), z_LQR(ii,:)', t_SNOPT, Lsol, z_SNOPT, u_SNOPT, R, polys, params);
    u_LQR(ii) = u;
end

%% Print Times output
% 

fprintf('Found the optimal trajectory in                      %f seconds.\n',params.tSNOPT);
fprintf('Solved the ODE using optimal open-loop control in    %f seconds.\n',params.tODE);
fprintf('Solved the Riccati ODE to find the LQR controller in %f seconds.\n',params.tRiccati');
fprintf('Solved the ODE using LQR closed-loop control in      %f seconds.\n',params.tLQR);


%% Plot the outputs for SNOPT and the plain ODE solution
% Still no control applied at this point.

close all

% Plot trajectory
figure, hold on
plot(z_SNOPT(:,1),z_SNOPT(:,2),'Linewidth',2)
    plot(z_OL(:,1),z_OL(:,2),'--','Linewidth',2)
        plot(z_LQR(:,1),z_LQR(:,2),':','Linewidth',2)
title('Trajectory','Interpreter','Latex')
xlabel('$x$ (m)','Interpreter','Latex')
ylabel('$y$ (m)','Interpreter','Latex')
legend({'SNOPT','OL','LQR'},'Interpreter','Latex')
axis equal

% Plot control effort
figure, hold on
plot(t_SNOPT,u_SNOPT,'Linewidth',2); 
plot(params.Tapp(:,1),params.Tapp(:,2),'--','Linewidth',2);
plot(t_LQR,u_LQR,':','Linewidth',2); 
title('Control','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Torque (Nm)','Interpreter','Latex')
legend({'$\tau_{SNOPT}$','$\tau_{ODE45}$'},'Interpreter','Latex')

% Plot cartesian positions
figure, hold on
plot(t_SNOPT,z_SNOPT(:,1),'Linewidth',2);  
    plot(t_OL, z_OL(:,1),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,1),':','Linewidth',2);
plot(t_SNOPT,z_SNOPT(:,2),'Linewidth',2);  
    plot(t_OL, z_OL(:,2),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,2),':','Linewidth',2);
title('Cartesian Positions','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Position (m)','Interpreter','Latex')
legend({'$x_{SNOPT}$','$x_{OL}$','$x_{LQR}$', ...
        '$y_{SNOPT}$','$y_{OL}$','$y_{LQR}$'},...
            'Interpreter','Latex','Location','northwest');
    
% Plot angular positions
figure, hold on
plot(t_SNOPT,z_SNOPT(:,3),'Linewidth',2);  
    plot(t_OL, z_OL(:,3),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,3),':','Linewidth',2);
plot(t_SNOPT,z_SNOPT(:,4),'Linewidth',2);  
    plot(t_OL, z_OL(:,4),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,4),':','Linewidth',2);
title('Angular Positions','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Angle (rad)','Interpreter','Latex')
legend({'$\theta_{SNOPT}$','$\theta_{OL}$','$\theta_{LQR}$', ...
        '$\phi_{SNOPT}$','$\phi_{OL}$','$\phi_{LQR}$'},...
            'Interpreter','Latex','Location','northeast');

% Plot cartesian velocities
figure, hold on
plot(t_SNOPT,z_SNOPT(:,5),'Linewidth',2);  
    plot(t_OL, z_OL(:,5),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,5),':','Linewidth',2);
plot(t_SNOPT,z_SNOPT(:,6),'Linewidth',2);  
    plot(t_OL, z_OL(:,6),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,6),':','Linewidth',2);
title('Cartesian Velocities','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Velocity (m/s)','Interpreter','Latex')
legend({'$\dot{x}_{SNOPT}$','$\dot{x}_{OL}$','$\dot{x}_{LQR}$', ...
        '$\dot{y}_{SNOPT}$','$\dot{y}_{OL}$','$\dot{y}_{LQR}$'},...
            'Interpreter','Latex','Location','northwest');

% Plot angular velocities
figure, hold on
plot(t_SNOPT,z_SNOPT(:,7),'Linewidth',2);  
    plot(t_OL, z_OL(:,7),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,7),':','Linewidth',2);
plot(t_SNOPT,z_SNOPT(:,8),'Linewidth',2);  
    plot(t_OL, z_OL(:,8),'--','Linewidth',2);
        plot(t_LQR, z_LQR(:,8),':','Linewidth',2);
title('Angular Velocities','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('$\omega$ (rad/s)','Interpreter','Latex')
legend({'$\dot{\theta}_{SNOPT}$','$\dot{\theta}_{OL}$','$\dot{\theta}_{LQR}$', ...
        '$\dot{\phi}_{SNOPT}$','$\dot{\phi}_{OL}$','$\dot{\phi}_{LQR}$'},...
            'Interpreter','Latex','Location','northeast');

function [dw] = dynamics(t, w, ts, Ls, z_states, z_controls, R, polys, params)

    [~,u] = findK(t, w, ts, Ls, z_states, z_controls, R, polys, params);
    dw = f(w,u,params);
    
%     waitbar(t/params.maxTime,params.waitBar,'Solving ODE with LQR')

end
 
