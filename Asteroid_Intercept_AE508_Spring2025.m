%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%

clear;
format longg;

% Calculate the state and costate differential equations
function Xdot = eom(t,X)

    r = X(1:3);
    v = X(4:6);
    m = X(7);
    %lam = X()

    rdot = v;

    % Xdot = [rdot; vdot; mdot];
end

function err = costFunction(t,X,T,c)
    % err = [];
end

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
% tf = 2462237.520833333; % A.D. 2029-Apr-11 00:30:00.0000 as Julian Date
% TODO: guess the time of flight that makes sense given the 
% satellite we want to hit the ateroid with is on a GEO orbit
tf = ??? % time of flight, in seconds
rf = [-1.096812308683544e+06, -9.292488001004120e+05, -4.114552797159857e+05]'; % the final position
% Velocity is free at final time
% vf = [4.206778300537007E+00, 3.801730069720992E+00, 1.639412680210433E+00]'; % the final velocity
% mass is free at final time
