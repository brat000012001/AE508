function err = max_velocity(lam0guess,t0,tf,x0,xf,T,c,rho,opts_ode,m0,m_a,mu)

    [~, X] = ode45(@eom, [t0 tf], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    % See my notes, page 120, boundary conditions
    PhiDot = terminal_cost(X(end,4:6),X(end,7),xf(4:6),m_a);
    err = [X(end,1:3)' - xf(1:3); X(end,11:14)' - PhiDot'];

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
