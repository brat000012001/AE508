function [p,f,g,h,k,L] = classical2mee(a,ecc,Om,incl,w,ta)
    %CLASSICAL2MEE Converts the classical orbital elements to
    % the Mixed Equinoctial Orbital Elements
    p = a*(1-ecc^2);
    f = ecc*cos(w+Om);
    g = ecc*sin(w+Om);
    h = tan(incl/2)*cos(Om);
    k = tan(incl/2)*sin(Om);
    L = Om + w + ta;
end