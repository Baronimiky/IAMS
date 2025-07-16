function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, om_f] = biellipticTransfer(ai, ei, af, ef, ra_t, om_i, om_f, P, mu)

% Bielliptic transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, om_f] = biellipticTransfer(ai, ei, af, ef, ra_t, om_i, om_f, P, mu)
%
% Input arguments:
% ai       [1x1]       Initial semi-major axis [km]
% ei       [1x1]       Initial eccentricity [-]
% af       [1x1]       Final semi-major axis [km]
% ef       [1x1]       Final eccentricity [-]
% ra_t     [1x1]       Apocenter of transfer orbit [km]
% om_i     [1x1]       Argument of pericenter (initial) [rad]
% om_f     [1x1]       Argument of pericenter (final) [rad]
% P        [1x1]       Starting point: 1=pericenter, 2=apocenter
% mu       [1x1]       Gravitational parameter [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1  [1x1]  1st maneuver impulse            [km/s]
% DeltaV2  [1x1]  2nd maneuver impulse            [km/s]
% DeltaV3  [1x1]  3nd maneuver impulse            [km/s]
% Deltat1  [1x1]  maneuver time                   [s]
% Deltat2  [1x1]  maneuver time                   [s]
% om_f     [1x1]  final pericenter anomaly        [rad]
%

% Calcola parametri orbitali
p1 = ai * (1 - ei^2);
p2 = af * (1 - ef^2);

% Rendi gli angoli positivi
om_f = PositiveAngles(om_f);
om_i = PositiveAngles(om_i);

% Calcola la differenza angolare
delta_omega = round(abs(om_i - om_f), 12);
isOpposite = (delta_omega == round(pi, 12));
isSame = (delta_omega == 0 || delta_omega == round(2*pi, 12));

% Determina r1 e r2
if (P == 2 && isOpposite)
    r1 = p1 / (1 - ei);  % apocentro iniziale
    r2 = p2 / (1 + ef);  % pericentro finale
end

if (P == 1 && isOpposite)
    r1 = p1 / (1 + ei);  % pericentro iniziale
    r2 = p2 / (1 - ef);  % apocentro finale
end

if (P == 2 && isSame)
    r1 = p1 / (1 - ei);  % apocentro iniziale
    r2 = p2 / (1 - ef);  % apocentro finale
end

if (P == 1 && isSame)
    r1 = p1 / (1 + ei);  % pericentro iniziale
    r2 = p2 / (1 + ef);  % pericentro finale
end

% Calcola semiassi maggiori delle orbite di trasferimento
at1 = (r1 + ra_t) / 2;
at2 = (ra_t + r2) / 2;

% Calcola i deltav
DeltaV1 = sqrt(2*mu*(1/r1 - 1/(2*at1))) - sqrt(2*mu*(1/r1 - 1/(2*ai)));
DeltaV2 = sqrt(2*mu*(1/ra_t - 1/(2*at2))) - sqrt(2*mu*(1/ra_t - 1/(2*at1)));
DeltaV3 = sqrt(2*mu*(1/r2 - 1/(2*af))) - sqrt(2*mu*(1/r2 - 1/(2*at2)));

% Calcola i tempi di trasferimento
Deltat1 = pi * sqrt(at1^3 / mu);
Deltat2 = pi * sqrt(at2^3 / mu);

% Determina l'argomento del pericentro finale
if ((isSame && P == 1) || (isOpposite && P == 2))
    om_f = om_i;  % nessuna rotazione
else
    om_f = PositiveAngles(om_i + pi);  % rotazione di 180°
end

end