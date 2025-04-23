function PhiDot = terminal_cost(v,m,va,ma)
%{
 Compute the partial derivatives of the terminal cost wrt state vector.
%}
    % The negative sign is intentional since we want to maximize
    % the momentum transfer
    denom = -sqrt(m^2*dot(v,v) + ma^2*dot(va,va) - 2*m*ma*dot(v,va));
    xdot = m*(m*v(1) - ma*va(1))/denom;
    ydot = m*(m*v(2) - ma*va(2))/denom;
    zdot = m*(m*v(3) - ma*va(3))/denom;
    mdot = (m*dot(v,v) - ma*dot(va,v))/denom;
    PhiDot = [xdot, ydot, zdot, mdot];
end