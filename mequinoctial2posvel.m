function [r_vec,v_vec] = mequinoctial2posvel(p,f,g,h,k,L,mu)
    %MEQUINOCTIAL2POSVEL Converts Mixed Equinoctial Elements to 
    % the position and velocity vectors in ECI/CBI

    s2 = 1 + h^2 + k^2;
    alpha2 = h^2 - k^2;
    w = 1 + f*cos(L) + g*sin(L);
    r = p/w;
    
    r_vec = [
            r/s2*(cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
            r/s2*(sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
            2*r/s2*(h*sin(L) - k*cos(L))
        ];
    v_vec = [
            -1/s2*sqrt(mu/p)*(sin(L) + alpha2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha2*g);
            -1/s2*sqrt(mu/p)*(-cos(L) + alpha2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha2*f);
            2/s2*sqrt(mu/p)*(h*cos(L) + k*sin(L) + f*h + g*k)
        ];
end