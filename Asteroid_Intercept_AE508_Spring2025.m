%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%

clearvars; close all; clc;
format longg;
addpath(".");

function write_to_file(filename,vec)
    fid = fopen(filename , 'a');
    fprintf(fid, [repmat('%.12f,\t', 1, size(vec,1)) '\n'], vec);
    fclose(fid);
end

function XDot = eom_aroid(t,X,mu)
    r = X(1:3);
    v = X(4:6);
    rdot = v;
    vdot = -mu/norm(r)^3*r;
    XDot = [rdot; vdot];
end

function XDot = eom_scraft(t,X,mu)
    XDot = eom_aroid(t,X,mu);
end

%
% Compute the partial derivative of the terminal cost wrt state vector
% The terminal cost is defined as diff
function PhiDot = terminal_cost(v,m,va,ma)
    denom = sqrt(m^2*dot(v,v) + ma^2*dot(va,va) - 2*m*ma*dot(v,va));
    xdot = -m*(m*v(1) - ma*va(1))/denom;
    ydot = -m*(m*v(2) - ma*va(2))/denom;
    zdot = -m*(m*v(3) - ma*va(3))/denom;
    mdot = -m*dot(v,v) - ma*dot(va,v);
    PhiDot = [xdot; ydot; zdot; mdot];
end

function err = minT(lam0guess,t0,tf,x0,xf,T,c,rho,opts_ode,m0,m_a,mu)
    tf_guess = lam0guess(8);
    [t, X] = ode45(@eom, [t0 tf_guess], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    % PhiDot = terminal_cost(X(end,4:6),X(end,7),xf(4:6),m_a);
    % See my notes, page 120, boundary conditions
    % err = [X(end,1:3)' - xf(1:3); X(end,11:14)' - PhiDot];
    H = hamiltonian(t,X,T,c,rho,mu);
    err = [X(end,1:3)' - xf(1:3); X(end,11:14)'; H(end) + 1];
end


state_values = init();

x0 = [state_values.r0;state_values.v0];
xf = [state_values.rf;state_values.vf];

rho = 1.0;
%
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);
%
% Solve Maximum momentum transfer optimal trajectory problem.
%
lam_guess = [1e-5*ones(7,1);state_values.tf];
iter = 0;
maxIter = 1000;
plotSet = plots();
while iter < maxIter
    write_to_file(state_values.logfile,lam_guess);
    iter = iter + 1;
    [p0,~,exitflag, output] = fsolve(@minT,...
        lam_guess, ...
        options, ...
        state_values.t0, ...
        state_values.tf, ...
        x0, ...
        xf, ...
        state_values.T, ...
        state_values.c, ...
        rho, ...
        opts_ode, ...
        state_values.m0, ...
        state_values.m_a, ...
        state_values.mu);
    if exitflag == 0
        fprintf("The number of iterations has been exceeded some threshold: %s\n", output.message);
        lam_guess = p0;
        continue;
    elseif exitflag < 0
        fprintf("The problem not solved: %s\n", output.message);
        break;
    else
        fprintf("The solution converged: %s\n", output.message);
        break;
    end
end
[t, X] = ode45(@eom, [state_values.t0 p0(8)], ...
    [x0; state_values.m0; p0(1:7)], ...
    opts_ode, ...
    state_values.T, ...
    state_values.c, ...
    rho, ...
    state_values.mu);

H = hamiltonian(t,X, ...
    state_values.T, ...
    state_values.c, ...
    rho, ...
    state_values.mu);
%
% Plot the Hamiltonian
%
plot(t, H, 'LineWidth', 1.5);
% Plot the states, costates, control and the optimal trajectory
plotSet.trajectory(t,X,state_values.r0,state_values.rf);

% Propagate the asteroid trajectory
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
[ta, Xa] = ode45(@eom_aroid, [0 86400*1.2], ...
    [state_values.rf; state_values.vf], ...
    opts_ode, ...
    state_values.mu);
[tsp, Xsp] = ode45(@eom_scraft, [0 86400], ...
    [state_values.r0; state_values.v_geo], ...
    opts_ode, ...
    state_values.mu);
hold on;
plot3(Xa(:,1),Xa(:,2),Xa(:,3),'LineWidth',1.5)
plot3(Xsp(:,1),Xsp(:,2),Xsp(:,3),'LineWidth',1.5)
hold off;

plotSet.states(t,X);
plotSet.costates(t,X);
plotSet.control(t,X);
