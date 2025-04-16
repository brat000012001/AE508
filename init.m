%
% Returns the initial and terminal state values (positions & velocities),
% thrust magnitude, the standard gravitational parameter, initial and  
% final times, the initial wet mass of the spacecraft, specific impulse,
%
function state = init()
    
    state = {};
    state.mu = 398600.44; % the Earth's gravitational parameter 
    
    % 11 Apr 2029 00:26:00.000 UTCG is approximately the date/time when asteroid Apophis 
    % enters the Earth's SOI.The SOI was calculated using the Hill Radius (~ 234 Earth radii)
    
    %
    % Initial condition
    % T = 1d, a = (T*sqrt(mu_earth)/(2*pi))^(2/3) = 42241.095610673336 km
    state.t0 = 0;
    state.r0 = [42241.095610673336, 0, 0]'; % km, GEO orbit
    state.v0 = [0, 3.0718591585665633, 0]'; % km/s the velocity of a GEO satellite
                                      % computed using the vis-viva eqn
    state.m0 = 500; % kg, the initial wet mass
    
    % Final condition
    % Asteroid Apophis position in Inertial reference frame centered at the
    % Earth center (ECI)
    state.rf = [-1.096812308683544e+06, -9.292488001004120e+05, -4.114552797159857e+05]'; % the final position
    % The asteroid velocity at the time of the intercept
    state.vf = [4.206778300537007E+00, 3.801730069720992E+00, 1.639412680210433E+00]'; % the final velocity
    % The mass of the hypothetical asteroid
    state.m_a = 2.2e8; % kg, the mass of 2024 YR4 asteroid
    
    % The time of flight that makes sense given the 
    % satellite we want to hit the ateroid with is on a GEO orbit
    state.tf = 3779250.236767169; % time of flight, in seconds
    
    state.opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
    state.options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
        'MaxIter',1e3,'TolFun',1e-12,'TolX',1e-14,...
        'UseParallel',false);
    
    state.T = 2.5035/1000; % kN, Thrust magnitude
    state.Isp = 4190; % seconds, specific impulse
    state.g0 = 9.8; % m/s^2 (Earth's surface gravity)
    state.c = state.Isp*state.g0/1000; % km/s Exhaust velocity
end