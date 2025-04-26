function err = max_dvast(lam0guess,t0,tf,x0,xf,T,c,rho,opts_ode,m0,m_ast,mu,beta)
    [~, X] = ode45(@eom, [t0 tf], [x0; m0; lam0guess(1:7)],opts_ode,T,c,rho,mu);
    % See my notes, page 120, boundary conditions
    PhiDot = terminal_cost(X(end,4:6)',X(end,7),xf(4:6));
    err = [X(end,1:3)' - xf(1:3); X(end,11:14)' - PhiDot'];

    function PhiDot = terminal_cost(v,m,v_ast)
    %{
     Compute the partial derivatives of the terminal cost wrt state vector.
    %}
        % The negative sign is intentional since we want to maximize
        % the momentum transfer
        dv = v - v_ast;
        DV2 = dv'*dv;
        denom = -(abs(beta)*m_ast*sqrt(DV2));
        xdot = (beta^2*m*dv(1))/denom;
        ydot = (beta^2*m*dv(2))/denom;
        zdot = (beta^2*m*dv(3))/denom;
        mdot = beta^2*DV2/denom;
        PhiDot = [xdot, ydot, zdot, mdot];
    end
end
