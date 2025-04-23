function XDot = two_body(t,X,mu)
%{
    Two-body dynamics
%}
    r = X(1:3);
    v = X(4:6);
    rdot = v;
    vdot = -mu/norm(r)^3*r;
    XDot = [rdot; vdot];
end