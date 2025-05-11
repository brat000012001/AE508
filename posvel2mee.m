function [p,f,g,h,k,L] = posvel2mee(r,v,cbi)
    %RV2MEE Convert the position and velocity vectors to mixed equinoctial elements
    % Arguments:
    % r - 3x1 position vector
    % v - 3x1 velocity vector
    % cbi - 3x3 central body inertial coordinate frame axes
    [a,ecc,Om,incl,w,ta] = posvel2classical(r,v,cbi);
    [p,f,g,h,k,L] = classical2mee(a,ecc,Om,incl,w,ta);
end