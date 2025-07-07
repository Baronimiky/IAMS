function [a,e_n,i,OM, om, th] = car2par(rr,vv,mu)

% Transformation from cartesian coordinates to Keplerian parameters
%
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% rr    [3x1]   position vector              [km]
% vv    [3x1]   velocity vector              [km/s]
% mu    [1x1]   gravitational parameter      [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% a     [1x1]   semi-major axis              [km]
% e     [1x1]   eccentricity                 [-]
% i     [1x1]   inclination                  [rad]
% OM    [1x1]   RAAN                         [rad]
% om    [1x1]   pericenter anomaly           [rad]
% th    [1x1]   true anomaly                 [rad]

r_n = norm(rr);
v_n = norm(vv);

E = 0.5*(v_n)^2 - mu/r_n;
a = (2/r_n -(v_n^2)/mu)^(-1);
h = cross(rr,vv);
h_n = norm(h);

e = cross(vv,h)/mu - rr/r_n;
e_n = norm(e);
i = acos (h(3)/h_n);
k = [0 0 1];

N = cross(k,h)/norm(cross(k,h));

if N(2) >=0
    OM = acos(N(1));
else
    OM = 2*pi - acos(N(1));
end

if e(3) >=0
    om = acos(dot(N,e)/e_n);
else
    om = 2*pi - acos(dot(N,e)/e_n);
end

vr = dot(vv,rr)/r_n;

if vr >= 0
    th = acos((dot(rr,e))/(r_n*e_n));
else
    th = acos((dot(rr,e))/(r_n*e_n));













end

