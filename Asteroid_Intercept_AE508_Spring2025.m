%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%
% troubleshooting:
% X(340169,4:6) NaN

clear; close all; clc;
format longg;
addpath(".");

% Calculate the state and costate differential equations
function Xdot = eom(t,X,T,c,rho,mu)

    r = X(1:3);
    v = X(4:6);
    m = X(7);
    lam_r = X(8:10);
    lam_v = X(11:13);
    lam_m = X(14);

    % S = c/m*norm(lam_v) + lam_m - 1; % switch function
    % delta = 0.5*(1 + tanh(S/rho));

    uhat = -lam_v/norm(lam_v);
    
    rdot = v;
    vdot = -mu/norm(r)^3*r + T/m*uhat;
    mdot = -T/c;

    x = r(1);
    y = r(2);
    z = r(3);

    scalar = -mu/((x^2+y^2+z^2)^(5/2));
    lam_r_dot = [
        scalar*lam_v'*[2*x^2-y^2-z^2;      3*x*y;            3*x*z];
        scalar*lam_v'*[3*x*y; -x^2 + 2*y^2 - z^2;            3*y*z];
        scalar*lam_v'*[3*x*z;              3*y*z;2*z^2 - x^2 - y^2]
        ];

    lam_v_dot = -lam_r;

    lam_m_dot = -T/m^2*(lam_v'*uhat);

    Xdot = [rdot; vdot; mdot; lam_r_dot; lam_v_dot; lam_m_dot];
end

function err = minT(lam0guess,t0,x0,xf,T,c,rho,opts_ode,m0,mu)
    tf = lam0guess(8);
    [t, X] = ode45(@eom, [t0 tf], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    H = hamiltonian(t,X,T,c,rho,mu);
    err = [X(end,1:3)' - xf(1:3); lam0guess(4:7); H(end) + 1];
end

function Hamiltonian = hamiltonian(t,X,T,c,rho,mu)
    Hamiltonian = zeros(length(t),1);
    rvec = X(:,1:3);
    vvec = X(:,4:6);
    mvec = X(:,7);
    rmag = vecnorm(rvec,2,2);
    lam_r = X(:,8:10);
    lam_v = X(:,11:13);
    lam_m = X(:, 14);
    uhat = -lam_v./vecnorm(lam_v,2,2);
    for idx = 1:length(t)
        Hamiltonian(idx) = lam_r(idx,:)*vvec(idx,:)' ...
            - lam_v(idx,:)*rvec(idx,:)'*mu/rmag(idx)^3 ...
            + lam_v(idx,:)*uhat(idx,:)'*T/mvec(idx) - lam_m(idx)*T/c; 
    end
end

mu_e = 398600.44; % the Earth's gravitational parameter 

% 11 Apr 2029 00:26:00.000 UTCG is approximately the date/time when asteroid Apophis 
% enters the Earth's SOI.The SOI was calculated using the Hill Radius (~ 234 Earth radii)

%
% Initial condition
% T = 1d, a = (T*sqrt(mu_earth)/(2*pi))^(2/3) = 42241.095610673336 km
t0 = 0;
r0 = [42241.095610673336, 0, 0]'; % km, GEO orbit
v0 = [0, 3.0718591585665633, 0]'; % km/s the velocity of a GEO satellite
                                  % computed using the vis-viva eqn
m0 = 500; % kg, the initial wet mass

% Final condition
% Asteroid Apophis position in Inertial reference frame centered at the
% Earth center (ECI)
rf = [-1.096812308683544e+06, -9.292488001004120e+05, -4.114552797159857e+05]'; % the final position
% Velocity is free at final time
vf = [4.206778300537007E+00, 3.801730069720992E+00, 1.639412680210433E+00]'; % the final velocity
% mass is free at final time

% The time of flight that makes sense given the 
% satellite we want to hit the ateroid with is on a GEO orbit
% t_min = 3335700.31527787
tf = 3336000.31527787; % time of flight, in seconds

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

T = 2.5035/1000; % kN, Thrust magnitude
Isp = 4190; % seconds, specific impulse
g0 = 9.8; % m/s^2 (Earth's surface gravity)
c = Isp*g0/1000; % km/s Exhaust velocity

% Solve the Lambert's equation analytically
[v1,v2] = lambert(r0,rf,tf,mu_e);

x0 = [r0;v1];
xf = [rf;0;0;0];

rho = 1.0;
%
lam_guess = [1e-5*ones(7,1);tf];
[p0,~] = fsolve(@minT,lam_guess,options,t0,x0,xf,T,c,rho,opts_ode,m0,mu_e);
[t, X] = ode45(@eom, [t0 p0(8)], [x0; m0; p0(1:7)],opts_ode,T,c,rho,mu_e);
H = hamiltonian(t,X,T,c,rho,mu_e);
plot(t,H);
plot_trajectory(t,X,r0,rf);
plot_states(t,X);
plot_costates(t,X);
plot_control(t,X);

%{
[v1,v2] = lambert(r0,rf,tf,mu_e);
lam_guess = [1e-5;1e-5;1e-5;1e-5;1e-5;1e-5;1e-5;tf];
[t, X] = ode45(@eom, [t0 tf], [[r0;v1]; m0; lam_guess(1:7)],opts_ode,T,c,rho,mu_e);
plot_trajectory(t,X,r0,rf);
%}

function plot_control(t,X)
    u = X(:,11:13)./vecnorm(X(:,11:13),2,2);

    figure;
    title('Control');
    subplot(3,1,1);
    plot(t,u(:,1),'LineWidth', 1.5);
    xlabel('t (sec)');
    ylabel('u_x');

    subplot(3,1,2);
    plot(t,u(:,2),'LineWidth', 1.5);
    xlabel('t (sec)');
    ylabel('u_y');
    
    subplot(3,1,3);
    plot(t,u(:,3),'LineWidth', 1.5);
    xlabel('t (sec)');
    ylabel('u_z');
end

function plot_trajectory(t, X, r0, rf)
    figure;
    % grid on; hold on;
    plot3(X(:,1),X(:,2),X(:,3));
    %scatter(0,0,0,'*r'); % The origin of the central body (the Earth)
    %scatter(r0(1),r0(2),r0(3),'*b');
    %plot3(rf(1),rf(2),rf(3),'*g');
    % hold off;
end

function plot_costates(t, X)
    figure;
    grid on; hold on;

    title('Costate \lambda_r');
    subplot(3,1,1);
    plot(t, X(:,8),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_r_x');

    subplot(3,1,2);
    plot(t, X(:,9),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_r_y');

    subplot(3,1,3);
    plot(t, X(:,10),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_r_z');
    hold off;
    
    % lambda_v
    figure;
    grid on; hold on;
    title ('Costate \lambda_v');
    subplot(3,1,1);
    plot(t, X(:,11),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_v_x');

    subplot(3,1,2);
    plot(t, X(:,12),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_v_y');

    subplot(3,1,3);
    plot(t, X(:,13),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_v_z');
    hold off;

    % Costate lambda_m
    figure;
    grid on; hold on;
    title('\lambda_m');
    plot(t, X(:,14),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('\lambda_m');
    hold off;
end

function plot_states(t,X)
    figure;
    grid on; hold on;

    title("Position")
    subplot(3,1,1);
    plot(t, X(:,1),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('x (km)');

    subplot(3,1,2);
    plot(t, X(:,2),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('y (km)');
    
    subplot(3,1,3);
    plot(t, X(:,3),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('z (km)');

    hold off;

    % Velocity
    figure;
    grid on; hold on;
    
    title("Velocity")
    subplot(3,1,1);
    plot(t, X(:,4),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('v_x (km/s)');

    subplot(3,1,2);
    plot(t, X(:,5),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('v_y (km/s)');
    
    subplot(3,1,3);
    plot(t, X(:,6),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('v_z (km/s)');
    
    hold off;

    % Mass 
    figure;
    grid on; hold on;

    title("Mass")
    plot(t, X(:,7),'LineWidth',1.5);
    xlabel('t (sec)');
    ylabel('m (kg)');
    
    hold off;
end