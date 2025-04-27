%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% The script performs a trade-off study by solving maximum momentum 
% trajectory problem for a satellite on an intercept trajectory with 
% some hypothetical  asteroid that is on a nar collision course to Earth
% by iterating over a range of initial spacecraft positions. 
% At each iteration the script solves the optional control trajectory
% problem by maximing the momentum transfer, computes the 
% transfer momentum at the final intercept position, and plots
% the results.
%
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%

clearvars; close all; clc;
format longg;
addpath(".");

%
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

filename = "reports/tradeoff_studies_results.txt";
ar = (0:10:360)';
times = flipud([4 5 6 7 8 9 10]');
thrusts = flipud([0.05 0.1 0.25 0.5 1.0 1.5 2]');

figure(1); hold on;
for timeidx = 1:length(times)
    time_of_flight = times(timeidx);
    for thrustidx = 1:length(thrusts)
        thrust_magnitude = thrusts(thrustidx);
        for idx = 1:length(ar)
        
            true_anomaly = deg2rad(ar(idx));
            state_values = initial_values(true_anomaly,time_of_flight,0.0, thrust_magnitude);
        
            x0 = [state_values.r0;state_values.v0];
            xf = [state_values.rf;state_values.vf];
        
            rho = 1.0;
        
            %
            % Solve Maximum momentum transfer optimal trajectory problem.
            %
            lam_guess = 1e-5*ones(7,1);
            [p0,~,exitflag, output] = fsolve(@max_momentum,...
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
        
            if exitflag > 0
        
                [~, X] = ode45(@eom, [state_values.t0 state_values.tf], ...
                    [x0; state_values.m0; p0(1:7)], ...
                    opts_ode, ...
                    state_values.T, ...
                    state_values.c, ...
                    rho, ...
                    state_values.mu);
                
                dv = sqrt( X(end,7)^2*norm(X(end,4:6))^2 ...
                    + state_values.m_a^2*norm(state_values.vf)^2 ...
                    - 2*X(end,7)*state_values.m_a*dot(X(end,4:6),state_values.vf));
        
                M = [dv, state_values.tf, true_anomaly, thrust_magnitude, x0', xf', p0'];
        
                if exist(filename, "file")
                    writematrix(M,filename,"WriteMode","append");
                else
                    writematrix(M,filename);
                end
                plot3(X(:,1),X(:,2),X(:,3));
            else
                M = [NaN, state_values.tf, true_anomaly, thrust_magnitude, x0', xf', p0'];
                if exist(filename, "file")
                    writematrix(M,filename,"WriteMode","append");
                else
                    writematrix(M,filename);
                end
            end
        end
    end
end
hold off;
%{
figure; hold on; grid on;
plot(ar, dv,'b*','LineWidth', 1.5);
xlabel('True Anomaly \nu (deg)');
ylabel('Transfer Momentum Magnitude');
hold off;
%}