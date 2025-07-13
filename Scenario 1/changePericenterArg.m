function [DeltaV, omf, thi, thf, DeltaT] = changePericenterArg (a, e, omi, omf, th_i, P, mu)
% Change of Pericenter Argument maneuver
%
% [DeltaV, thetai, thetaf] = changePericenterArg (a, e, omi, omf, mu)
%
% ---------------------------------------------------------------------
%
% Input arguments:
% a         [1x1]       semi-major axis                 [km]
% e         [1x1]       eccentricity                    [-]
% omi       [1x1]       initial pericenter anomaly      [rad]
% omf       [1x1]       final pericenter anomaly        [rad]
% mu        [1x1]       gravitational parameter         [km^3/s^2]
%
% ---------------------------------------------------------------------
% Output arguments:
% DeltaV    [1x1]       maneuver impulse                [km/s]
% omf       [1x1]       final pericenter anomaly        [km/s]
% thi       [1x1]       initial true anomaly            [rad]
% thf       [1x1]       final true anomaly              [rad]
%


%  con P=1: 
%    se a è il punto di manovra più vicino (vicino al perigeo) e b è il punto di manovra più
%    lontano (vicino all'apogeo), allora se thi è dopo a e prima di b la
%    manovra avviene in b, altrimenti avviene in a
% con P=2: 
%    prende l'altro punto (quello più lontano rispetto al senso di
%    movimento del corpo lungo l'orbita)


omf=PositiveAngles(omf);
omi=PositiveAngles(omi);

dw = omf-omi;


if (dw/2)<th_i && th_i<=(pi+dw/2)
    thf=pi-dw/2;
    thi=pi+(dw/2);
else
    thi = (dw)/2;
    thf = 2 * pi - (dw)/2;
end

if P==2
    thf=thf+pi;
    thi=thi+pi;
end

p = a * (1 - e^2);

DeltaV = 2 * sqrt(mu / p) * e * sin((dw)/2);
DeltaT = TOF(a,e,th_i,thi,mu);

end