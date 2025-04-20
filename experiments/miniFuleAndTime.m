%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Author: Peter Nalyvayko (petern4@illinois.edu)
% 
% Trajectory optimization 
%

clearvars; close all; clc;
format longg;
addpath("..");

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


function err = minT(lam0guess,t0,x0,xf,T,c,rho,opts_ode,m0,mu)
    tf_guess = lam0guess(8);
    [t, X] = ode45(@eom, [t0 tf_guess], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    H = hamiltonian(t,X,T,c,rho,mu);
    % See my notes, page 101, the last term H(end) + 1
    % v and m are free parameters, so their costates must be 0
    err = [X(end,1:3)' - xf(1:3); X(end,11:14)'; H(end) + 1];
end

state_values = init();

x0 = [state_values.r0;state_values.v0];
xf = [state_values.rf;state_values.vf];

rho = 1.0;
%
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',2e3,...
    'MaxIter',2e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

lam_guess = [1e-5*ones(7,1); state_values.tf];
[p0,~,exitflag, output] = fsolve(@minT, ...
    lam_guess, ...
    options, ...
    state_values.t0, ...
    x0, ...
    xf, ...
    state_values.T, ...
    state_values.c, ...
    rho, ...
    opts_ode, ...
    state_values.m0, ...
    state_values.mu);

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
plotSet = plots();
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
