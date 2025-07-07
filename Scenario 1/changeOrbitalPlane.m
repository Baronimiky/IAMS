function [DeltaV, omf, theta]=changeOrbitalPlane(a,e,i_i,OMi,omi,i_f,OMf,thi,P,mu)

% Change of Plane maneuver
%
% [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a       [1x1]    semi-major axis                   [km]
% e       [1x1]    eccentricity                      [-]
% i_i     [1x1]    initial inclination               [rad]
% OMi     [1x1]    initial RAAN                      [rad]
% omi     [1x1]    initial pericenter anomaly        [rad]
% i_f     [1x1]    final inclination                 [rad]
% OMf     [1x1]    final RAAN                        [rad]
% mu      [1x1]    gravitational parameter           [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV  [1x1]    maneuver impulse                  [km/s]
% omf     [1x1]    final pericenter anomaly


if nargin<10
    mu=398600.433;
end

alpha=acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(OMf-OMi));
ui=atan2(sin(OMf-OMi)/sin(alpha)*sin(i_f), (-cos(i_f)+cos(alpha)*cos(i_i))/sin(alpha)/sin(i_i));
uf=atan2(sin(OMf-OMi)/sin(alpha)*sin(i_i), (cos(i_i)-cos(alpha)*cos(i_f))/sin(alpha)/sin(i_f));


if((OMf-OMi)>0 && (i_f-i_i)>0)
    thi = ui - omi;
    omf = uf - thi;
elseif((OMf-OMi)>0 && (i_f-i_i)<0)
    thi = 2*pi - ui - omi;
    omf = 2 * pi - uf - thi;
elseif((OMf-OMi)<0 && (i_f-i_i)>0)
    alpha=acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(abs(OMf-OMi)));
    ui=atan2(sin(abs(OMf-OMi))/sin(alpha)*sin(i_f), (-cos(i_f)+cos(alpha)*cos(i_i))/sin(alpha)/sin(i_i));
    uf=atan2(sin(abs(OMf-OMi))/sin(alpha)*sin(i_i), (cos(i_i)-cos(alpha)*cos(i_f))/sin(alpha)/sin(i_f));
    thi = 2*pi - ui - omi;
    omf = 2 * pi - uf - thi;
elseif((OMf-OMi)<0 && (i_f-i_i)<0)
    alpha=acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(abs(OMf-OMi)));
    ui=atan2(sin(abs(OMf-OMi))/sin(alpha)*sin(i_f), (-cos(i_f)+cos(alpha)*cos(i_i))/sin(alpha)/sin(i_i));
    uf=atan2(sin(abs(OMf-OMi))/sin(alpha)*sin(i_i), (cos(i_i)-cos(alpha)*cos(i_f))/sin(alpha)/sin(i_f));
    thi = ui - omi;
    omf = uf - thi;
end
theta = thi;

if cos(theta)>0
    theta = theta + pi;
end

DeltaV=2*sqrt(mu/(a*(1-e^2)))*(1+e*cos(theta))*sin(alpha/2);

end