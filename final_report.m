%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% The script generates the plots for the final project report.
%
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%
clearvars; close all; clc;
format longg;
addpath(".");

rho = 1.0;

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

results = readmatrix("reports/tradeoff_studies_results_final.txt");
solutions = results(~isnan(results(:,1)),:);
% dv (1)
% tf (2)
% true_anomaly (3)
% thrust_magnitude (4)
% x0' (5:10)
% xf' (11:16)
% p0' (17:23)

default_state = initial_values();

times = unique(solutions(:,2));
max_by_tof_thrust = [];
max_velocity = [];
for timeidx = 1:size(times,1)
    figure(timeidx); hold on;
    by_time = solutions(solutions(:,2) == times(timeidx),:);
    thrusts = unique(by_time(:,4));
    rows = idivide(size(thrusts,1),int16(2));
    if mod(size(thrusts,1),2) > 0
        rows = rows + 1;
    end

    rows = double(rows);
    for thrustidx = 1:length(thrusts)
        if rows > 1 || length(thrusts) > 1
            subplot(rows,2,thrustidx); hold on;
        end
        title(sprintf("Thrust: %.1f N, ToF: %.2f days", thrusts(thrustidx), times(timeidx)/86400));
        xlabel('x (km)');
        ylabel('y (km)');
        zlabel('z (km)');
        plot3(0,0,0,'b*','LineWidth',3);
        plot3(default_state.rf(1),default_state.rf(2),default_state.rf(3),'r*','LineWidth',3);
        
        by_thrust = by_time(by_time(:,4) == thrusts(thrustidx),:);
         [~,maxidx] = max(by_thrust(:,1));
         max_by_tof_thrust(size(max_by_tof_thrust,1)+1,:) = by_thrust(maxidx,:);
        for rowidx = 1:size(by_thrust,1)
            x = by_thrust(rowidx,:);
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
    
            plot3(X(1,1),X(1,2),X(1,3),'g+','LineWidth',1);

            if rowidx == maxidx
                max_velocity(size(max_velocity,1)+1,:) = [by_thrust(maxidx,1), by_thrust(maxidx,4), X(end,4:6)];
                plot3(X(:,1),X(:,2),X(:,3),'r--','LineWidth',1.5);
            else
                plot3(X(:,1),X(:,2),X(:,3));
            end
        end
    end
    hold off;
end
%{
x = solutions(:,2)/86400;
y = solutions(:,4);
z = solutions(:,1);
xi = unique(x);
yi = unique(y);
[X,Y] = meshgrid(xi,yi);
Z = reshape(z, size(X));
figure;
surf(X,Y,Z);
%}