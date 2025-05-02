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

spacecraft_velocity_at_impact = [];
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
    spacecraft_velocity_at_impact(size(spacecraft_velocity_at_impact,1)+1,:) =...
        [idx,v_tf',m_tf,x(4),x(2)];
end

asteroid_masses = linspace(2.2e6,2.2e8,50);
asteroid_deltavs = [];

for massidx = 1:length(asteroid_masses)
    m_ast = asteroid_masses(massidx);
    for idx = 1:size(spacecraft_velocity_at_impact,1)
        x = spacecraft_velocity_at_impact(idx,:);
        dv = deltavast(x(2:4)',x(5), default_values.vf,m_ast,1.0);
        % asteroid_deltavs(i,:) => 
        % 1: index of the original solution, 
        % 2-4: new asteroid velocity vector, 3x1
        % 5: new asteroid velocity magnitude
        % 6: spacecraft thrust magnitude
        % 7: time of flight
        % 8: asteroid mass
        asteroid_deltavs(size(asteroid_deltavs,1)+1,:) = [x(1),...
            (default_values.vf + dv)',...
            norm(dv),...
            x(6),...
            x(7),...
            m_ast];
    end
end

% Propagate forward the asteroid orbit using the position
% and velocity of the asteroid at the time of impact in ECI
[~, Xast] = ode45(@two_body, [0 86400], [default_values.rf;default_values.vf], ...
    opts_ode, default_values.mu);
ast_to_earth_distances = vecnorm(Xast(:,1:3),2,2);
ast_to_earth_min_distance = min(ast_to_earth_distances);

% Propagate forward the asteroid orbit using the estimated
% change in the asteroid's velocity vector due to impact, 
% compute the closest approach of the asteroid on a new orbit
% to Earth
final_results = [];
for idx = 1:size(asteroid_deltavs,1)
    x = asteroid_deltavs(idx,:);
    vast_new = x(2:4)';
    [~, Xastnew] = ode45(@two_body, [0 86400], [default_values.rf;vast_new], ...
        opts_ode, default_values.mu);
    Xastdist = vecnorm(Xastnew(:,1:3),2,2);
    closest_approach = min(Xastdist);
    if closest_approach > ast_to_earth_min_distance
        final_results(size(final_results,1)+1,:) = [x, closest_approach];
    end
end
%
% Sort by closest distance
sorted = sortrows(final_results,9,'descend');
best_solutions = sorted(1:min(4,size(sorted,1)),:);

figure; hold on;
% Plot the first four solutions 
for bestidx = 1:size(best_solutions,1)
    row = best_solutions(bestidx,:);
    x = solutions(row(1),:);
    [~,~,X] = propagate(x, opts_ode);
    plot3(X(:,1),X(:,2),X(:,3)); hold on;
end

msg1 = sprintf("dv=%g km/s, T=%.1f N, t_f=%.1f days, m_ast=%.3f", ...
    best_solutions(1,5),best_solutions(1,6),best_solutions(1,7)/86400,best_solutions(1,8));
msg2 = sprintf("dv=%g km/s, T=%.1f N, t_f=%.1f days, m_ast=%.3f", ...
    best_solutions(2,5),best_solutions(2,6),best_solutions(2,7)/86400,best_solutions(2,8));
msg3 = sprintf("dv=%g km/s, T=%.1f N, t_f=%.1f days, m_ast=%.3f", ...
    best_solutions(3,5),best_solutions(3,6),best_solutions(3,7)/86400,best_solutions(3,8));
msg4 = sprintf("dv=%g km/s, T=%.1f N, t_f=%.1f days, m_ast=%.3f", ...
    best_solutions(4,5),best_solutions(4,6),best_solutions(4,7)/86400,best_solutions(4,8));

% Plot the Earth and asteroid positions
plot3(0,0,0, 'g*','LineWidth',3);
plot3(default_values.rf(1),default_values.rf(2),default_values.rf(3),'r*','LineWidth',3);

title('Top 4 best solutions');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
legend(msg1,msg2,msg3,msg4,'Earth center', 'Asteroid @ t_f', 'Location','best');
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

function dv = deltavast(v,m,v_ast,m_ast,beta)
    dv = beta*m/m_ast*(v - v_ast);
end