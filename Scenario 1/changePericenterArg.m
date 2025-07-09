function [DeltaV, thi, thf] = changePericenterArg (a, e, omi, omf, mu)
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
% thi       [2x1]       initial true anomalies          [rad]
% thf       [2x1]       final true anomalies            [rad]
%

omf=PositiveAngles(omf);
omi=PositiveAngles(omi);

dw = omf-omi;
theta1i = (dw)/2;
theta1f = 2 * pi - (dw)/2;

theta2i = pi + (dw)/2;
theta2f = pi - (dw)/2;

thi = [theta1i; theta2i];
thf = [theta1f; theta2f];

p = a * (1 - e^2);

DeltaV = 2 * sqrt(mu / p) * e * sin((dw)/2);

end