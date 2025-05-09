%
% AE508 Optimal Space Trajectories, Spring 2025

% Example: Spin out trajectory

clear; close all; clc;
format longg;

% Calculate the state and costate differential equations
function Xdot = eom(t,X,T,c,mu)

    r = X(1:3);
    v = X(4:6);
    m = X(7);

    uhat = v/norm(v);
    
    rdot = v;
    vdot = -mu/norm(r)^3*r + T/m*uhat;
    mdot = -T/c;
    Xdot = [rdot; vdot; mdot];
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
% vf = [4.206778300537007E+00, 3.801730069720992E+00, 1.639412680210433E+00]'; % the final velocity
% mass is free at final time

% tf = 2462237.520833333; % A.D. 2029-Apr-11 00:30:00.0000 as Julian Date
theta = acos((r0/norm(r0))'*(rf/norm(rf)));
c = sqrt(r0'*r0 + rf'*rf - 2*norm(r0)*norm(rf)*cos(theta));
s = (norm(r0) + norm(rf) + c)/2;
% Minimum time of flight
t_p = sqrt(2)/(3*sqrt(mu_e))*(sqrt(s^3) - sign(sin(theta))*sqrt((s - c)^3));

% The time of flight that makes sense given the 
% satellite we want to hit the ateroid with is on a GEO orbit
tf = t_p; % time of flight, in seconds

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
    'UseParallel',false);

T = 1.5935/1000; % kN, Thrust magnitude
Isp = 4190; % seconds, specific impulse
g0 = 9.8; % m/s^2 (Earth's surface gravity)
c = Isp*g0/1000; % km/s Exhaust velocity

x0 = [r0;v0];
xf = [rf;0;0;0];

[t, X] = ode45(@eom, [t0 tf], [x0; m0],opts_ode,T,c,mu_e);
plot_trajectory(t,X,r0,rf);

function plot_trajectory(t, X, r0, rf)
    figure;
    plot3(X(:,1),X(:,2),X(:,3));
end

