%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% Author: Peter Nalyvayko (petern4@illinois.edu)
% Source: Orbital Mechanics for Engineering Students 
function [v1,v2] = lambert(p1, p2, t, mu, string)
%{
The function solves the Lambert's Problem.

p1, p2   - initial and final position vectors (km)
v1, v2   - initial and final velocity vectors (km/s)
mu       - gravitational parameter (km^3/s^2)
string   - 'pro'   if the trajectory is prograde
           'retro' if the trajectory is retrograde
%}

    r1 = norm(p1);
    r2 = norm(p2);
    
    c12 = cross(p1,p2);
    theta = acos(dot(p1,p2)/r1/r2);

    % Determine whether orbit is prograde or retrograde
    if nargin < 5 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))
        string = 'pro';
        fprintf("\n Prograde trajectory assumed.\n");
    end

    if strcmp(string,'pro')
        if c12(3) <= 0
           theta = 2*pi - theta;
        end
    elseif strcmp(string,'retro')
        if c12(3) >= 0
            theta = 2*pi - theta;
        end
    end

    % Equation 5.35
    A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

    % Determine approx where F(z,t) changes sign, and use
    % that value of z as the starting value for Equation 5.45
    z = -100;
    while F(z,t) < 0
        z = z + 0.1;
    end

    % Setan error tolerance and limit on the number of iterations
    tol = 1.e-8;
    nmax = 5000;

    % Iterate on Equation 5.45 until z is determined to within the error
    % tol
    ratio = 1;
    n = 0;
    while (abs(ratio) > tol) & (n < nmax)
        n = n + 2;
        ratio = F(z,t)/dFdz(z);
        z = z - ratio;
    end

    % Report if the maximum number of iterations is exceeded
    if n >= nmax
        fprintf('\n\n ** Number of iterations exeeded %g \n\n', nmax);
    end

    % Equation 5.46a
    f = 1 - y(z)/r1;

    % Equation 5.46b
    g = A*sqrt(y(z)/mu);

    % Equation 5.46d
    gdot = 1 - y(z)/r2;

    % Equation 5.28
    v1 = 1/g*(p2 - f*p1);

    % Equation 5.29
    v2 = 1/g*(gdot*p2 - p1);

    return;
    
    % Equation 5.38
    function dum = y(z)
        dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    end

    % Equation 5.40
    function dum = F(z,t)
        dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
    end

    % Equation 5.43
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z))...
                + 3*S(z)*2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z))...
                + A*sqrt(C(z)/y(z)));
        end
    end

    % Stumpff functions
    function dum = C(z)
        dum = stumpC(z);
    end

    function dum = S(z)
        dum = stumpS(z);
    end

    function s = stumpS(z)
        % Evaluate the Stumpff function S(z) according to
        % Equation 3.52
        if z > 0
            s = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
        elseif z < 0
            s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
        else
            s = 1/6;
        end
    end

    function c = stumpC(z)
        % Evaluates the Stumpff function C(z) according to Equation 3.53
        if z > 0
            c = (1 - cos(sqrt(z)))/z;
        elseif z < 0
            c = (cosh(sqrt(-z)) - 1)/(-z);
        else
            c = 1/2;
        end
    end

end % end lambert

