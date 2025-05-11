function [a,ecc,Om,incl,w,ta] = mee2classical(p,f,g,h,k,L)
    %MEE2CLASSICAL Converts the Modified Equinoctial Orbital Elements
    % to the Classical Orbital Elements
    a = p/(1 - f^2 - g^2);
    ecc = sqrt(f^2 + g^2);
    incl = atan2(2*sqrt(h^2+k^2),1 - h^2 - k^2);
    w = atan2(g*h - f*k, f*h + g*k);
    Om = atan2(k,h);
    ta = L - atan(g/f);
end