function [r_vec,v_vec] = classical2posvel(a,ecc,lan,incl,aop,ta, mu)
%{
    Orbital elements to position/velocity vectors using
    Direction Cosine Matrix approach
%}
    % Specific angular momentum: h^2/mu = a(1-e^2)
    h = sqrt(mu*a*(1-ecc^2));
    % orbit position magnitude
    r = h^2/(mu*(1+ecc*cos(ta)));
    % position and velocity vectors in perifocal frame
    r_pf = r*[cos(ta); sin(ta); 0];
    v_pf = mu/h*[-sin(ta); ecc+cos(ta); 0];
    dcm = DirectionCosineMatrix();
    r_vec = dcm*r_pf;
    v_vec = dcm*v_pf;

    function C = DirectionCosineMatrix()
      C1 = [cos(lan) sin(lan) 0; -sin(lan) cos(lan) 0; 0 0 1];
      C2 = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
      C3 = [cos(aop) sin(aop) 0; -sin(aop) cos(aop) 0; 0 0 1];
      C = transpose(C3*C2*C1);
    end
  
end
