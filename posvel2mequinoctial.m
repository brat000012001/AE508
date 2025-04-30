function [p,f,g,h,k,L] = posvel2mequinoctial(r,v,cbi)
    %RV2MEE Convert the position and velocity vectors to mixed equinoctial elements
    % Arguments:
    % r - 3x1 position vector
    % v - 3x1 velocity vector
    % cbi - 3x3 central body inertial coordinate frame axes
    [a,ecc,Om,incl,w,ta] = posvel2classical(r,v,cbi);
    p = a*(1-ecc^2);
    f = ecc*cos(w+Om);
    g = ecc*sin(w+Om);
    h = tan(incl/2)*cos(Om);
    k = tan(incl/2)*sin(Om);
    L = Om + w + ta;
end