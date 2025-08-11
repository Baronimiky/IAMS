function [dv,omf,thf,dt]=changeOrbitalPlane(ai,ei,Ii,OMi,omi,If,OMf,thi,P,mu)

% Change of Plane maneuver
%
% [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, thi, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai       [1x1]    semi-major axis                   [km]
% ei       [1x1]    eccentricity                      [-]
% Ii       [1x1]    initial inclination               [rad]
% OMi      [1x1]    initial RAAN                      [rad]
% omi      [1x1]    initial pericenter anomaly        [rad]
% If       [1x1]    final inclination                 [rad]
% OMf      [1x1]    final RAAN                        [rad]
% thi      [1x1]    initial point ~= manouver point   [rad]
% P - Maneuver point selection:
%       1: Perform maneuver at farthest point from current position
%       2: Perform maneuver at closest point to current position
% mu       [1x1]    gravitational parameter           [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV   [1x1]    maneuver impulse                  [km/s]
% omf      [1x1]    final pericenter anomaly          [rad]
% theta    [1x1]    true anomaly at manouver          [rad]
% dt       [1x1]    time from starting point to manouver point    [s]

Ii=PositiveAngles(Ii);
omi=PositiveAngles(omi);
OMi=PositiveAngles(OMi);
If=PositiveAngles(If);
OMf=PositiveAngles(OMf);

if nargin==9
    mu=398600;
end

alpha=acos(cos(Ii)*cos(If)+sin(Ii)*sin(If)*cos(OMf-OMi));


if((OMf-OMi)>=0 && (If-Ii)>=0)
    ui=acos((cos(alpha)*cos(Ii)-cos(If))/(sin(alpha)*sin(Ii)));
    uf=acos((cos(Ii)-cos(alpha)*cos(If))/(sin(alpha)*sin(If)));
    if (OMf-OMi)>pi
        ui=-ui;
        uf=-uf;
    end
end
if((OMf-OMi)>=0 && (If-Ii)<0) 
    ui=-acos((cos(If)-cos(alpha)*cos(Ii))/(sin(alpha)*sin(Ii)));
    uf=-acos((cos(alpha)*cos(If)-cos(Ii))/(sin(alpha)*sin(If)));
    if (OMf-OMi)>pi
        uf=-uf;
        ui=-ui;
    end
end
if((OMf-OMi)<0 && (If-Ii)>=0)
    ui=-acos((cos(alpha)*cos(Ii)-cos(If))/(sin(alpha)*sin(Ii)));
    uf=-acos((cos(Ii)-cos(alpha)*cos(If))/(sin(alpha)*sin(If)));
    if (OMf-OMi)<-pi
       ui=-ui;
       uf=-uf;
    end  
end
if((OMf-OMi)<0 && (If-Ii)<0) 
    ui=acos((cos(If)-cos(alpha)*cos(Ii))/(sin(alpha)*sin(Ii)));
    uf=acos((cos(alpha)*cos(If)-cos(Ii))/(sin(alpha)*sin(If)));
    if (OMf-OMi)<-pi
        uf=-uf;
        ui=-ui;
    end
end    

if(ui>=omi)
    thf=ui-omi;
else
    thf=2*pi+(ui-omi);
end

if(uf>=thf)
    omf=uf-thf;
else
    omf=2*pi+(uf-thf);
end

if (P==2)
    thf=thf+pi;
end

dv=2*sqrt(mu/(ai*(1-ei^2)))*(1+ei*cos(thf))*sin(alpha/2);

dt=TOF(ai,ei,PositiveAngles(thi),PositiveAngles(thf),mu);