function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)
% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = bitangentTransfer(ai, ei, af, ef, ra_t, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai        [1x1]   initial semi-major axis       [km]
% ei        [1x1]   initial eccentricity          [-]
% af        [1x1]   final semi-major axis         [km]
% ef        [1x1]   final eccentricity            [-]
% ra_t      [1x1]   transfer orbits apocenter distance [km]
% mu        [1x1]   gravitational parameter       [km^3/s^2]
%
% Output arguments:
% DeltaV1   [1x1]   1st maneuver impulse          [km/s]
% DeltaV2   [1x1]   2nd maneuver impulse          [km/s]
% DeltaV3   [1x1]   3rd maneuver impulse          [km/s]
% Deltat1   [1x1]   maneuver time 1               [s]
% Deltat2   [1x1]   maneuver time 2               [s]

rp_t1 = ai*(1-ei);
rp_t2 = af*(1-ef);

a_t1 = (rp_t1+ra_t)/2;
a_t2 = (rp_t2+ra_t)/2;

DeltaV1 = sqrt(mu)*sqrt((2/rp_t1)-(1/a_t1)) - sqrt(mu)*sqrt((2/rp_t1)-(1/ai));
DeltaV2 = sqrt(mu)*sqrt((2/ra_t)-(1/a_t2)) - sqrt(mu)*sqrt((2/ra_t)-(1/a_t1));
DeltaV3 = sqrt(mu)*sqrt((2/rp_t2)-(1/af)) - sqrt(mu)*sqrt((2/rp_t2)-(1/a_t2));

Deltat1 = 2*pi * sqrt((a_t1^3)/mu);
Deltat2 = 2*pi * sqrt((a_t2^3)/mu);


end