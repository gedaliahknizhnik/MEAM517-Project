function [z,F,INFO] = optimizeChaplygin(x_0, x_f, ind1, ind2, nx, nu, N, dt, params)
% OPTIMIZECHAPLYGIN(nx, nu, N, dt, params) uses SNOPT to solve the direct
%       collocation problem concerning the Chaplygin Beanie. Modified from
%       the HW 5 code creat by Mathew Halm.
% ------------------------------------------------------------------------
% INPUTS:
%   nx - the size of the state vector x
%   nu - the size of the control vector u
%   N  - the number of segments
%   dt - the time between segments (constant)
%   params - a structure of problem parameters (masses, inertias, lengths,
%   etc.)
% OUTPUTS
%   z - the final state vector (including all xs and us)
%   F - the final constraint values
%   INFO - the output code from SNOPT

%% Define problem settings:


%% Define Constraints
% Aeq*z = beq

Aeq = zeros(2*nx, N * (nx + nu));
beq = zeros(2*nx, 1);

% Add constraints to Aeq, beq to enforce 
% starting at x_0 and ending at x_f
x_0_inds = 1:nx;
xf_inds = x_0_inds + (N - 1) * (nx + nu);

% Place identity matrices in the x_0 block and the x_N block of Aeq.
% A1 = zeros(nx);
% A2 = zeros(nx);
% for ii = 1:nx
%     A1(ii,ii) = ind1(ii);
%     A2(ii,ii) = ind2(ii);
% end

% Aeq(x_0_inds,x_0_inds) = A1;
% Aeq(x_0_inds+nx,xf_inds) = A2;

%Only constrict position
% f_con = eye(nx);
% f_con(nx/2+1:end,nx/2+1:end) = zeros(nx/2);
% Aeq(x_0_inds+nx,xf_inds) = f_con;

    
% Let the first nx elements of beq be x_0 and the last nx elements be
% x_f.
% beq(x_0_inds)    = x_0;
% beq(x_0_inds+nx) = x_f;

%% Add Constraints on States

M = Inf;

lb = -inf(N * (nx + nu),1);
ub =  inf(N * (nx + nu),1);

% Add constraint on initial position
for i=1:nx
    if (ind1(i) == 1)
        lb(i) = x_0(i);
        ub(i) = x_0(i);
    end
    
    if (ind2(i) == 1)
        index = size(lb,1) - (nx+nu);
        lb(index + i) = x_f(i);
        ub(index + i) = x_f(i);
    end
end

% Add bounding box constraints u_1 \in [-M,M], u_2 \in [-M,M]
for i=1:N
  u_i_inds = (1:nu) + nx * i + nu * (i - 1);
  lb(u_i_inds) = -M;
  ub(u_i_inds) = M;
end
  
% Make initial guess for z
z0 = zeros(N * (nx + nu), 1);

for i=1:N
  x_i_inds = (1:nx) + (nx + nu) * (i - 1);
  u_i_inds = (1:nu) + nx * i + nu * (i - 1);
  z0(x_i_inds) = x_0 + ((i-1)/(N-1))*(x_f - x_0); 
  %z0(u_i_inds) = M;
end

%% Solve the problem using SNOPT

snscreen on;
snprint('chaplygin.out');
snsummary('chaplygin_summary.out')

snset('Verify')

[z,zlow,zupp,zmul,zstate, ...
   Flow,Fupp,Fmul,Fstate, ...
 ObjAdd,ObjRow,A,iAfun,jAvar,iGfun,jGvar,iii,jjj] = setOptData(z0,lb,ub);
       
[z,F,INFO,zmul,Fmul] = snopt(z, zlow, zupp, zmul, zstate, ...
		   Flow, Fupp, Fmul, Fstate, ...
		   @usrfun, ObjAdd, ObjRow, A, iAfun, jAvar, iGfun, jGvar);

snprint off; % Closes the file and empties the print buffer
snend;

%% Suplemental Functions

    % SNOPT takes in F = [g;h], where g is the cost function and h is the
    % vector of all the constraints. In the current iteration I am not
    % passing any gradients.
    
    function [F,dF] = usrfun(z) 
        [g, dG] = trajectory_cost(z, N, nx, nu, dt);
        [c,ceq,dC,dCeq] = all_constraints(z, N, nx, nu, dt, params);       
        F = [g; ceq];
    
        % Add the non-zero elements of the cost gradient.
        dGnonZero = zeros(N,1);
        for ii=1:N
            u_i_inds = (1:nu) + nx * ii + nu * (ii - 1);
            dGnonZero(ii) = dG(u_i_inds);
        end
        
        % Add the non-zero elements of the constraint gradient        
        s = zeros(size(iii,1),1);
        for ii = 1:size(iii,1)
            s(ii) = dCeq(iii(ii),jjj(ii));
        end

        % Combine all the gradients into a single vector
        dF = [dGnonZero;s];
    end
% 
    % We also have to set up all of the other variables that SNOPT requires
    function [z,zlow,zupp,zmul,zstate, ...
            Flow,Fupp,Fmul,Fstate, ...
            ObjAdd,ObjRow,A,iAfun,jAvar,iGfun,jGvar,iii,jjj] = setOptData(z0,lb,ub)

        A = []; iAfun=[]; jAvar=[];
        
        ObjRow = 1;
        ObjAdd = 0;

        z      = z0;
        zlow   = lb;
        zupp   = ub;
        zmul   = []; zstate = [];

        [g, dG] = trajectory_cost(z0, N, nx, nu, dt);
        [iii,jjj] = genSparseStructure(nx,nu,N);
        
        length = 1 + nx*(N-1);

        Flow    = zeros(length,1);
        Flow(1) = -Inf;
        Fupp    = zeros(size(Flow));
        Fupp(1) =  Inf; 
        Fmul   = []; Fstate = [];
                
        iGfun = zeros(N+size(iii,1),1);
        jGvar = zeros(size(iGfun,1),1);
        
        % Encode the gradient of the cost (first row of F) with respect to
        % each of the control inputs.
        for ii=1:N
            u_i_inds = (1:nu) + nx * ii + nu * (ii - 1);
            iGfun(ii) = 1;
            jGvar(ii) = u_i_inds;
        end
        
        % Encode the non-zero elements of the gradient of the constraints
        iGfun(N+1:end) = iii+1;
        jGvar(N+1:end) = jjj;

    end

end