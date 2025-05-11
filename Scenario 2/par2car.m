function [rr,vv] = par2car(a,e,i,OM,om,th,mu)

% Trasformation from Keplerian parameters
%
% [rr,vv] = par2car(a,e,i,OM,om,th,mu)
%
% --------------------------------------------------------
% Input arguments:
% a         [1x1]  semi-major axis




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

