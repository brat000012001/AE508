function [a,ecc,Om,incl,w,ta] = posvel2classical(r,v, cbi)
    %POSVEL2CLASSICAL Converts position/velocity vectors to Keplerian
    %elements
    %   Takes the position and velocity vectors and converts them to the
    %   six classical elements
    h = cross(r,v);
    ecc_vec = cross(v,h)/mu - r/norm(r);
    ecc = norm(ecc_vec);
    % a(1-ecc^2)=h^/mu
    a = dot(h,h)/(mu*(1-ecc^2));
    %
    % Inclination
    %
    e1 = cbi(:,1);
    e2 = cbi(:,2);
    e3 = cbi(:,3);
    n = cross(e3,h);% Line of nodes
    nxe3 = cross(n,e3);
    incl = atan2(dot(h,nxe3/norm(nxe3)),dot(h,e3)); 
    % Compute the longitude of ascending node
    Om = atan2(dot(n,e2),dot(n,e1)); % Omega
    % Argument of periapsis
    hxn = cross(h/norm(h),n/norm(n));
    w = atan2(dot(ecc_vec,hxn/norm(hxn)),dot(ecc_vec,n/norm(n)));
    % True anomaly
    hxe = cross(h,ecc_vec);
    ta = atan2(dot(r,hxe/norm(hxe)),dot(r,ecc_vec/norm(ecc_vec)));

    if ta < 0
      ta = ta + 2*pi;
    end

    w = mod(w, 2*pi);
    Om = mod(Om, 2*pi);
end