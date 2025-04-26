
function state = initial_values(trueAnomaly, timeOfFlight, incl, thrust_magnitude)
% Returns the initial and terminal state values (positions & velocities),
% thrust magnitude, the standard gravitational parameter, initial and  
% final times, the initial wet mass of the spacecraft, specific impulse,
%
% Arguments:

% - trueAnomaly:     the true anomaly orbital parameter for the spacecraft on GEO
%                    orbit. Specified in radians
% - timeOfFlight:    desired fixed time of flight. Specified in days.
% - incl:            departure orbit inclination (radians). Default is 0 rad
% - thrust_magnitude: Thrust Magnitude in N
%
    addpath(".");
    state.mu = 398600.44; % the Earth's gravitational parameter 
    %
    % Initial condition
    state.t0 = 0;
    state.v_geo = [0, 3.0718591585665633, 0]'; % km/s the velocity of a GEO satellite
                                      % computed using Vis-Viva Eqn

    R = 42241.095610673336; % the radius of a GEO orbit

    if nargin < 1
        trueAnomaly = deg2rad(10);
    end

    if nargin < 2
        timeOfFlight = 10; % by default, the time of flight is 10 days
    end

    if nargin < 3
        incl = 0;
    end

    if nargin < 4
        thrust_magnitude = 2.5;
    end

    [r0,v0] = classical2posvel(R, 0.0, 0.0, incl, 0.0, trueAnomaly, state.mu);
    state.r0 = r0; % km/s initial velocity
    state.v0 = sqrt(2*state.mu/R)*v0/norm(v0); % km/s initial velocity
    % The position of asteroid Apophis at 2462239.715277778, A.D. 2029-Apr-13 05:10:00.0000
    state.rf = [-2.956932341462374E+05, -2.043495023251738E+05, -9.847039118095844E+04]';
    state.vf = [4.292596028662857E+00,  3.878022383303825E+00,  1.677451136476388E+00]';
    state.tf = timeOfFlight * 86400; % optimal time of flight 5.29552915031233 days
    state.m0 = 500; % kg, the initial wet mass
    % The time of flight that makes sense given the 
    % satellite we want to hit the ateroid with is on a GEO orbit
    state.Isp = 4190; % seconds, specific impulse
    state.g0 = 9.8; % m/s^2 (Earth's surface gravity)
    state.c = state.Isp*state.g0/1000; % km/s Exhaust velocity
    state.T = thrust_magnitude/1000; % kN, Thrust magnitude
    % The mass of the hypothetical asteroid
    state.m_a = 2.2e8; % kg, the mass of 2024 YR4 asteroid
end