function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, om_f] = biellipticTransfer(ai, ei, af, ef, ra_t, om_i, mu)
% Bielliptic transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, om_f] = biellipticTransfer(ai, ei, af, ef, ra_t, mu, om_i)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai        [1x1]   Initial semi-major axis           [km]
% ei        [1x1]   Initial eccentricity              [-]
% af        [1x1]   Final semi-major axis             [km]
% ef        [1x1]   Final eccentricity                [-]
% ra_t      [1x1]   Apocenter of transfer orbit       [km]
% mu        [1x1]   Gravitational parameter           [km^3/s^2]
% om_i      [1x1]   Argument of pericenter (initial)  [rad]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1   [1x1]   1st maneuver impulse              [km/s]
% DeltaV2   [1x1]   2nd maneuver impulse              [km/s]
% DeltaV3   [1x1]   3rd maneuver impulse              [km/s]
% Deltat1   [1x1]   Time of first transfer orbit      [s]
% Deltat2   [1x1]   Time of second transfer orbit     [s]
% om_f      [1x1]   Argument of pericenter (final)    [rad]


    rp_t1 = ai * (1 - ei);
    rp_t2 = af * (1 - ef);

    a_t1 = (rp_t1 + ra_t) / 2;
    a_t2 = (rp_t2 + ra_t) / 2;

    DeltaV1 = sqrt(mu) * (sqrt((2 / rp_t1) - (1 / a_t1)) - sqrt((2 / rp_t1) - (1 / ai)));
    DeltaV2 = sqrt(mu) * (sqrt((2 / ra_t) - (1 / a_t2)) - sqrt((2 / ra_t) - (1 / a_t1)));
    DeltaV3 = sqrt(mu) * (sqrt((2 / rp_t2) - (1 / af)) - sqrt((2 / rp_t2) - (1 / a_t2)));

    Deltat1 = pi * sqrt(a_t1^3 / mu);  % half-period of first transfer
    Deltat2 = pi * sqrt(a_t2^3 / mu);  % half-period of second transfer

    % Determine maneuver type and compute om_f
    tol = 1e-6;
    if abs(rp_t1 - rp_t2) < tol
        % pp case
        om_f = om_i;
    elseif abs((2 * ra_t - rp_t1 - rp_t2)) < tol
        % aa case
        om_f = om_i;
    else
        % pa or ap case
        om_f = PositiveAngles(om_i + pi);
    end

end