function [a,e,i,OM, om, th] = car2par(rr,vv,mu)
% Trasformation from cartesian coordinates to Keplerian parameters
%[a, e, i, OM, om, th] = car2par (rr, vv, mu)
% Input arguments :
% rr[3x1] position vector   [ km ]
% VV [3x1] velocity vector  [ km/s]
% mu [1x1] gravitational parameter[ km^3/s^2]
%
%
% Output arguments:
% a [1x1] semi-major axis     [km]
% e [1x1] eccentricity        [-]
% i [1x1] inclination         [rad]
% OM[1x1] RAAN                [rad]
% om [1x1] pericenter anomaly [ rad]
% th [1x1] true anomaly       [rad]

% ---calculate norms of the state vectors
r_norm = norm(rr);
v_norm = norm(vv);

% ---radial velocity and position WRT the pericenter
vr = dot(rr,vv)./r_norm;

% ---find the semi-major axis from the energy formula
epsilon = (1/2).*(v_norm.^2) - mu./r_norm; 
a = - mu./(2 * epsilon);

% ---specific angular momentum
h = cross(rr, vv);
h_norm = norm(h);

% ---inclination in radians
i = acos(h(3)./h_norm);             % 0 < i < pi

% ---nodes line
K = [0, 0, 1];
N = cross(K,h);
N_norm = norm(N);

% taking care of the singularity
% the nodes line will coincide with the Vernal Equinox line
if i == 0
    OM = 0;
end

% ---longitude of ascending node
if N(2) >= 0
    OM = acos(N(1)/N_norm);
else 
    OM = 2*pi - acos (N(1)/N_norm); 
end 

% ---eccentricity
e_vect = (1./mu).*((v_norm.^2 -mu./r_norm).*rr - (vr*r_norm).*vv); 
e = norm(e_vect);

% ---argument of the periapsis
if (e_vect(3) >= 0)
    om = acos(dot(N,e_vect)./(N_norm.*e)); 
else 
    om = 2*pi - acos(dot(N,e_vect)./(N_norm.*e)); 
end 

% taking care of the singularity
% the nodes line will coincide with the eccentricity vector
if e == 0
    om = 0;
end

 % ---argument of the periapsis
if vr >= 0
    th = acos(dot(e_vect,rr)./(e.*r_norm)); 
else 
    th = 2*pi - acos(dot(e_vect,rr)./(e.*r_norm));
end  

end











