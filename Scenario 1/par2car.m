function [rr,vv] = par2car(a,e,i,OM,om,th,mu)

% Transformation from Keplerian parameters to cartesian coordinates
%
% [rr, vv] = par2car(a, e, i, OM, om, th, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a     [1x1]   semi-major axis              [km]
% e     [1x1]   eccentricity                 [-]
% i     [1x1]   inclination                  [rad]
% OM    [1x1]   RAAN                         [rad]
% om    [1x1]   pericenter anomaly           [rad]
% th    [1x1]   true anomaly                 [rad]
% mu    [1x1]   gravitational parameter      [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% rr    [3x1]   position vector              [km]
% vv    [3x1]   velocity vector              [km/s]


p = a*(1-e^2);

r_n = p/(1+e*cos(th));

r_temp = r_n* [cos(th); sin(th);0];
v_temp = sqrt(mu/p) * [-sin(th); e+cos(th); 0];

R_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T_eci_pf = R_om*R_i*R_OM;

rr = T_eci_pf' * r_temp;
vv = T_eci_pf' * v_temp;

end