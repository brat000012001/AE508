function err = max_momentum_mee(lam0guess,t0,tf,x0,xf,T,c,rho,opts_ode,m0,m_a,mu)

    [~, X] = ode45(@eom_mee, [t0 tf], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    % MEE to position/velocity
    h = X(end,1);
    k = X(end,2);
    p = X(end,3);
    f = X(end,4);
    g = X(end,5);
    L = X(end,6);
    m = X(end,7);
    [r,v] = mee2posvel(p,f,g,h,k,L,mu);
    % See my notes, page 120, boundary conditions
    PhiDot = terminal_cost(v,m,xf(4:6),m_a);
    err = [abs(r - xf(1:3)) < 1e-4; X(end,8:14)' - PhiDot'];

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
end
