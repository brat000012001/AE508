%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%

%
% Calculates the state and costate differential equations
%
function Xdot = eom(t,X,T,c,rho,mu)

    r = X(1:3);
    v = X(4:6);
    m = X(7);
    lam_r = X(8:10);
    lam_v = X(11:13);
    lam_m = X(14);
    R3 = norm(r)^3;
    R5 = norm(r)^5;

    % S = c/m*norm(lam_v) + lam_m - 1; % switch function
    % delta = 0.5*(1 + tanh(S/rho));

    u = -lam_v/norm(lam_v);
    
    rdot = v;
    vdot = -mu/R3*r + T/m*u;
    mdot = -T/c;

    x = r(1);
    y = r(2);
    z = r(3);

    g11 = 3*mu*x^2/R5 - mu/R3;
    g12 = 3*mu*x*y/R5;
    g13 = 3*mu*x*z/R5;
    g21 = 3*mu*y*x/R5;
    g22 = 3*mu*y^2/R5 - mu/R3;
    g23 = 3*mu*y*z/R5;
    g31 = 3*mu*z*x/R5;
    g32 = 3*mu*z*y/R5;
    g33 = 3*mu*z^2/R5 - mu/R3;
    G = [g11 g12 g13; g21 g22 g23; g31 g32 g33];

    % Page 136 (Lecture 9 Optimal CSI Trajectory), 
    % Page 245 (Lecture 16. Cooperative Rendezvous Problem), my notes
    lam_r_dot = transpose(-lam_v'*G);
    lam_v_dot = -lam_r;
    lam_m_dot = T/m^2*(lam_v'*u);

    Xdot = [rdot; vdot; mdot; lam_r_dot; lam_v_dot; lam_m_dot];
end
