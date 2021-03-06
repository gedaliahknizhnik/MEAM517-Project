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

x_0 = [0,0,0,0,0,0,0,0];
x_f = [5,0,0,0,0,0,0,0];

ind1 = [1,1,1,1,1,1,1,1];
ind2 = [1,1,1,0,0,0,1,1];

M = Inf;

% Setup simulation parameters.
params.m     = 2;
params.B     = 1e0;
params.C     = 1e0;
params.a     = 0.05;
params.Tresb = 0;
params.Tresc = 0;
params.Fx    = 0;
params.Fy    = 0;
params.Fxw   = 0;
params.Fyw   = 0;
params.maxTime = T;
params.k1 = 0.0;
params.nx = nx;
params.nu = nu;

%% Solve the Optimal Control Problem

tic
[z,F,INFO] = optimizeChaplygin(x_0,x_f,ind1,ind2,M,nx,nu,N,dt,params);
toc

num = size(z,1)/(nx+nu);

z_states = zeros(num,nx);
z_controls = zeros(num,nu);

nx_indices = unwrap(1:nx);
nu_indices = nx + unwrap(1:nu);

for ii = 1:num
    x_indices = (ii-1)*(nx+nu) + nx_indices;
    u_indices = (ii-1)*(nx+nu) + nu_indices;
    z_states(ii,:)   = z(x_indices)';
    z_controls(ii,:) = z(u_indices)';
end

t = zeros(size(z_states,1),1);
for ii=1:size(t)
    t(ii) = (ii-1)*dt;
end

max(F)

%% Simulate ODE with these inputs 

Tapp = zeros(size(z_controls,1),2);

for ii = 1:size(Tapp,1)
    Tapp(ii,:) = [(ii-1)*dt, z_controls(ii)];
end

params.Tapp = Tapp;

W = waitbar(0,'Solving ODE'); params.waitBar = W;
[t_cont,z_cont] = ode45(@chaplyginSleigh,[0,T],x_0,[],params);
close(W)


%% Plot the outputs

close all
figure, hold on
plot(z_states(:,1),z_states(:,2),'--','Linewidth',2)
plot(z_cont(:,1),z_cont(:,2),'Linewidth',2)
title('Trajectory','Interpreter','Latex')
xlabel('$x$ (m)','Interpreter','Latex')
ylabel('$y$ (m)','Interpreter','Latex')
legend({'SNOPT','ODE45'},'Interpreter','Latex')
axis equal

figure, hold on
plot(t,z_states(:,1),'--','Linewidth',2);  plot(t_cont, z_cont(:,1),'Linewidth',2);
plot(t,z_states(:,2),'--','Linewidth',2);  plot(t_cont, z_cont(:,2),'Linewidth',2);
plot(t,z_states(:,3),'--','Linewidth',2);  plot(t_cont, z_cont(:,3),'Linewidth',2); 
plot(t,z_states(:,4),'--','Linewidth',2);  plot(t_cont, z_cont(:,4),'Linewidth',2);
title('Positions','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value ($m$,$rad$)','Interpreter','Latex')
legend({'$x_{SNOPT}$','$x_{ODE45}$','$y_{SNOPT}$','$y_{ODE45}$',...
            '$\theta_{SNOPT}$','$\theta_{ODE45}$','$\phi_{SNOPT}$',...
            '$\phi_{ODE45}$'},'Interpreter','Latex')

figure, hold on
plot(t,z_states(:,5),'--','Linewidth',2);  plot(t_cont, z_cont(:,5),'Linewidth',2);
plot(t,z_states(:,6),'--','Linewidth',2);  plot(t_cont, z_cont(:,6),'Linewidth',2);
plot(t,z_states(:,7),'--','Linewidth',2);  plot(t_cont, z_cont(:,7),'Linewidth',2);
plot(t,z_states(:,8),'--','Linewidth',2);  plot(t_cont, z_cont(:,8),'Linewidth',2);
title('Velocities','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value ($\frac{m}{s}$,$\frac{rad}{s}$)','Interpreter','Latex')
legend({'$\dot{x}_{SNOPT}$','$\dot{x}_{ODE45}$','$\dot{y}_{SNOPT}$',...
            '$\dot{y}_{ODE45}$','$\dot{\theta}_{SNOPT}$',...
            '$\dot{\theta}_{ODE45}$','$\dot{\phi}_{SNOPT}$',...
            '$\dot{\phi}_{ODE45}$'},'Interpreter','Latex')

figure, hold on
plot(t,z_controls,'--','Linewidth',2); plot(params.Tapp(:,1),params.Tapp(:,2),'Linewidth',2);
title('Control','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value (Nm)','Interpreter','Latex')
legend({'$\tau_{SNOPT}$','$\tau_{ODE45}$'},'Interpreter','Latex')

%% LQR


Q = 1*eye(nx);
Q(nx/2+1:end,nx/2+1:end) = 0.1*eye(nx/2);

R = 100*eye(nu);
F = zeros(size(Q));
tol = 1e-6;  % Accuracy of ricatti propagation

polys = getPolys(z_states,z_controls,nx,dt,params);

Soln = trajectoryLqr(t,z_states,z_controls,@linSys,Q,R,F,tol,polys,params);


ts = t;

Lf = reshape(0.01*eye(nx),nx^2,1);
myfun = @(t,L) Ldot(t, L, Q, R, ts, z_states, z_controls, polys, params);
[ts1,Ls1] = ode45(myfun,[T,0], Lf);

size(ts1)
size(Ls1)

return


ts = zeros(size(Soln,2),1);
Ks = zeros(size(Soln,2),nx);

for ii=1:size(Soln,2)
    ts(ii) = Soln(ii).t;
    Ks(ii,:) = Soln(ii).K;
end

W = waitbar(0,'Solving ODE with LQR'); params.waitBar = W;
[t_lqr,z_lqr] = ode45(@dynamics,[0,T],x_0,[],ts,Ks,z_states,z_controls, polys, params);
close(W)

u_lqr = zeros(size(t_lqr));

for ii = 1:size(t_lqr)
    u_star = interp1(ts,z_controls,t_lqr(ii));
    w_star = interp1(ts,z_states,t_lqr(ii))';
    KNow = interp1(ts,Ks,t_lqr(ii));
    
    u_lqr(ii) = u_star + KNow*(w_star-z_lqr(ii,:)');
end

%% Plot the outputs

figure, hold on
plot(z_states(:,1),z_states(:,2),'--','Linewidth',2)
plot(z_lqr(:,1),z_lqr(:,2),'Linewidth',2)
title('Trajectory','Interpreter','Latex')
xlabel('$x$ (m)','Interpreter','Latex')
ylabel('$y$ (m)','Interpreter','Latex')
legend({'SNOPT','ODE45/LQR'},'Interpreter','Latex')
axis equal

figure, hold on
plot(t,z_states(:,1),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,1),'Linewidth',2);
plot(t,z_states(:,2),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,2),'Linewidth',2);
plot(t,z_states(:,3),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,3),'Linewidth',2); 
plot(t,z_states(:,4),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,4),'Linewidth',2);
title('Positions','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value ($m$,$rad$)','Interpreter','Latex')
legend({'$x_{SNOPT}$','$x_{LQR}$','$y_{SNOPT}$','$y_{LQR}$',...
            '$\theta_{SNOPT}$','$\theta_{LQR}$','$\phi_{SNOPT}$',...
            '$\phi_{LQR}$'},'Interpreter','Latex')

figure, hold on
plot(t,z_states(:,5),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,5),'Linewidth',2);
plot(t,z_states(:,6),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,6),'Linewidth',2);
plot(t,z_states(:,7),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,7),'Linewidth',2);
plot(t,z_states(:,8),'--','Linewidth',2);  plot(t_lqr, z_lqr(:,8),'Linewidth',2);
title('Velocities','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value ($\frac{m}{s}$,$\frac{rad}{s}$)','Interpreter','Latex')
legend({'$\dot{x}_{SNOPT}$','$\dot{x}_{LQR}$','$\dot{y}_{SNOPT}$',...
            '$\dot{y}_{LQR}$','$\dot{\theta}_{SNOPT}$',...
            '$\dot{\theta}_{LQR}$','$\dot{\phi}_{SNOPT}$',...
            '$\dot{\phi}_{LQR}$'},'Interpreter','Latex')
        
figure, hold on
plot(t,z_controls,'--','Linewidth',2); plot(t_lqr,u_lqr,'Linewidth',2);
title('Control','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('Value (Nm)','Interpreter','Latex')
legend({'$\tau_{SNOPT}$','$\tau_{LQR}$'},'Interpreter','Latex')




function [dw] = dynamics(t, w, ts, Ks, z_states, z_controls, polys, params)
    
%    [w_star,u_star] = getDes(t,ts,z_states,z_controls,polys);

     u_star = interp1(ts,z_controls,t);
     w_star = interp1(ts,z_states,t)';

    KNow = interp1(ts,Ks,t);
    
    u = u_star + KNow*(w_star-w);
    dw = f(w,u,params);
    
    waitbar(t/params.maxTime,params.waitBar,'Solving ODE with LQR')

end
 
