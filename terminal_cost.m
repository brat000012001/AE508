function PhiDot = terminal_cost(v,m,va,ma)
%{
 % Compute the partial derivative of the terminal cost wrt state vector
 % The terminal cost is defined as diff
%}
%
    denom = -sqrt(m^2*dot(v,v) + ma^2*dot(va,va) - 2*m*ma*dot(v,va));
    xdot = m*(m*v(1) - ma*va(1))/denom;
    ydot = m*(m*v(2) - ma*va(2))/denom;
    zdot = m*(m*v(3) - ma*va(3))/denom;
    mdot = (m*dot(v,v) - ma*dot(va,v))/denom;
    PhiDot = [xdot, ydot, zdot, mdot];
%{
    V3 = norm(v)^3;
    V2 = norm(v)^2;
    VA = norm(va);
    VA2 = norm(va)^2;
    denom = -V3*VA*sqrt(V2 + VA2 - 2*dot(v,va)/V2/VA2);
    PhiDot(1) = v(1)*V3*VA - v(3)^2*va(1) - v(2)^2*va(1) + v(1)*v(2)*va(2) + v(1)*v(3)*va(3);
    PhiDot(2) = v(2)*V3*VA - v(3)^2*va(2) - v(1)^2*va(2) + v(1)*v(2)*va(1) + v(2)*v(3)*va(3);
    PhiDot(3) = v(3)*V3*VA - v(2)^2*va(3) - v(1)^2*va(3) + v(1)*v(3)*va(1) + v(2)*v(3)*va(2);
    PhiDot(4) = 0;
    PhiDot = PhiDot./denom;
%}
end