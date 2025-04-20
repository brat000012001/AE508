%
% Returns the initial and terminal state values (positions & velocities),
% thrust magnitude, the standard gravitational parameter, initial and  
% final times, the initial wet mass of the spacecraft, specific impulse,
%
function state = init(experimentid)
%{
    experimentid   - 0: (default) pick the location of asteroid Apophis 
                        at approx 373527 km from Earth. This is the default
                     1: pick the location of asteroid Apophis at the edge
                        of the Earth's SOI (calculated as Hills' radius at 234 @
                        Earth radii)
%}
    
    if nargin < 1 || experimentid == 0
        state = defaults();
    elseif experimentid == 1
        state = experimentSOI();
    end

    state.mu = 398600.44; % the Earth's gravitational parameter 
    %
    % Initial condition
    % T = 1d, a = (T*sqrt(mu_earth)/(2*pi))^(2/3) = 42241.095610673336 km
    state.t0 = 0;
    state.v_geo = [0, 3.0718591585665633, 0]'; % km/s the velocity of a GEO satellite
                                      % computed using Vis-Viva Eqn
    state.r0 = [42241.095610673336, 0, 0]'; % km, GEO orbit
    state.v0 = [0, 4.34426488374484, 0]'; % km/s initial velocity
    state.Isp = 4190; % seconds, specific impulse
    state.g0 = 9.8; % m/s^2 (Earth's surface gravity)
    state.c = state.Isp*state.g0/1000; % km/s Exhaust velocity
    % The mass of the hypothetical asteroid
    state.m_a = 2.2e8; % kg, the mass of 2024 YR4 asteroid
    state.m0 = 500; % kg, the initial wet mass
    state.T = 2.5/1000; % kN, Thrust magnitude

    function state = experimentSOI()
        state = {};
        % The position of asteroid Apophis at the edge of the Earth's SOI
        % at approx A.D. 2029-Apr-11 00:26:00.00. The coordinates are relative
        % the ecliptic of J2000
        % 11 Apr 2029 00:26:00.000 UTCG is approximately the date/time when asteroid Apophis 
        % enters the Earth's SOI.The SOI was calculated using the Hill Radius (~ 234 Earth radii)
        state.rf = [-1.096812308683544e+06, -9.292488001004120e+05, -4.114552797159857e+05]';
        state.vf = [4.206778300537007E+00, 3.801730069720992E+00, 1.639412680210433E+00]';
        % The time of flight that makes sense given the 
        % satellite we want to hit the ateroid with is on a GEO orbit
        state.tf = 3779250.236767169; % time of flight, in seconds

        state.logfile = "initial_guess_1.txt";
    end
    
    function state = defaults()
        state = {};
        % The position of asteroid Apophis at 2462239.715277778, A.D. 2029-Apr-13 05:10:00.0000
        state.rf = [-2.956932341462374E+05, -2.043495023251738E+05, -9.847039118095844E+04]';
        state.vf = [4.292596028662857E+00,  3.878022383303825E+00,  1.677451136476388E+00]';
        state.m0 = 500; % kg, the initial wet mass
        % The time of flight that makes sense given the 
        % satellite we want to hit the ateroid with is on a GEO orbit
        state.tf = 9.68335648148148 * 86400; % time of flight, in seconds

        state.logfile = "initial_guess_default.txt";
    end
end