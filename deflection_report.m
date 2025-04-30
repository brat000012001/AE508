%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% The script loads the optimal trajectory solutions generated
% by tradeoff_studies.m script, and plots top ten solutions
% with the highest chance of deflecting the asteroid.
%
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%
clearvars; close all; clc;
format longg;
addpath(".");

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

% Load the results
results = readmatrix("reports/tradeoff_studies_results_final.txt");
% Filter the results to keep only the solutions that converged
solutions = results(~isnan(results(:,1)),:);

default_values = initial_values();

asteroid_masses = [2.2e2 2.2e3 2.2e4 2.2e5 2.2e6 2.2e7 2.2e8];
asteroid_deltavs = [];
for idx = 1:size(solutions,1)
    x = solutions(idx,:);
    % dv (1)
    % tf (2)
    % true_anomaly (3)
    % thrust_magnitude (4)
    % x0' (5:10)
    % xf' (11:16)
    % p0' (17:23)
    p0 = x(17:23)';
    x0 = x(5:10)';
    [v_tf,m_tf,~] = propagate(x, opts_ode);
    [dv,dvmag] = deltavast(v_tf,m_tf, default_values.vf, default_values.m_a,1.0);
    asteroid_deltavs(size(asteroid_deltavs,1)+1,:) = [idx,dvmag,x(4)];
end

sorted = sortrows(asteroid_deltavs,2,'descend');
% Plot the first ten solutions with the highest change in the asteroid
% velocity after impact
figure; hold on;
for bestidx = 1:3
    row = asteroid_deltavs(bestidx,:);
    x = solutions(row(1),:);
    [~,~,X] = propagate(x, opts_ode);
    plot3(X(:,1),X(:,2),X(:,3));
end
msg1 = sprintf("dv=%.6f km/s, T=%.2f",asteroid_deltavs(1,2),asteroid_deltavs(1,3));
msg2 = sprintf("dv=%.6f km/s, T=%.2f",asteroid_deltavs(2,2),asteroid_deltavs(2,3));
msg3 = sprintf("dv=%.6f km/s, T=%.2f",asteroid_deltavs(3,2),asteroid_deltavs(3,3));
legend(msg1,msg2,msg3);
hold off;

function [v_tf,m_tf,X] = propagate(x, opts_ode)
    rho = 1.0;
    % dv (1)
    % tf (2)
    % true_anomaly (3)
    % thrust_magnitude (4)
    % x0' (5:10)
    % xf' (11:16)
    % p0' (17:23)
    p0 = x(17:23)';
    x0 = x(5:10)';
    state_values = initial_values(x(3),x(2)/86400,0.0,x(4));
    [~, X] = ode45(@eom, [state_values.t0 state_values.tf], ...
        [x0; state_values.m0; p0(1:7)], ...
        opts_ode, ...
        state_values.T, ...
        state_values.c, ...
        rho, ...
        state_values.mu);
    v_tf = X(end,4:6)';
    m_tf = X(end, 7);
end

function [dv,dvmag] = deltavast(v,m,v_ast,m_ast,beta)
    dv = beta*m/(m + m_ast)*(v - v_ast);
    dvmag = norm(dv);
end