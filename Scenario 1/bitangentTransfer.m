function [DeltaV1, DeltaV2, Deltat, om_f] = bitangentTransfer(ai, ei, af, ef, om_i, type, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, Deltat, om_f] = bitangentTransfer(ai, ei, af, ef, om_i, type, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai        [1x1]  initial semi-major axis          [km]
% ei        [1x1]  initial eccentricity             [-]
% af        [1x1]  final semi-major axis            [km]
% ef        [1x1]  final eccentricity               [-]
% om_i      [1x1]  argument of pericenter           [rad]
% type      [char] maneuver type   ('pa', 'ap', 'pp', 'aa')
% mu        [1x1]  gravitational parameter          [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1   [1x1]  1st maneuver impulse            [km/s]
% DeltaV2   [1x1]  2nd maneuver impulse            [km/s]
% Deltat    [1x1]  maneuver time                   [s]
% omf       [1x1]  final pericenter anomaly        [rad]
%

    Deltat = 0;

    % calcolo i raggi di pericentro e apocentro iniziali e finali
    rp_i = ai * (1 - ei);
    ra_i = ai * (1 + ei);
    rp_f = af * (1 - ef);
    ra_f = af * (1 + ef);

    switch type
        case 'pa'  % Pericentro a Apocentro
            rpT = rp_i;
            raT = ra_f;
        case 'ap'  % Apocentro a Pericentro
            rpT = rp_f;
            raT = ra_i;
        case 'pp'  % Pericentro a Pericentro
            rpT = rp_i;
            raT = rp_f;
        case 'aa'  % Apocentro a Apocentro
            rpT = ra_i;
            raT = ra_f;
        otherwise
            error('Scelta non valida. Possibilità: pa, ap, pp, aa.');
    end
    
    % Semi-asse maggiore dell'orbita di trasferimento
    aT = (rpT + raT) / 2;
    
    % Velocità
    vp_i = sqrt(mu * (2 / rp_i - 1 / ai));
    va_i = sqrt(mu * (2 / ra_i - 1 / ai));
    vp_f = sqrt(mu * (2 / rp_f - 1 / af));
    va_f = sqrt(mu * (2 / ra_f - 1 / af));
    vpT = sqrt(mu * (2 / rpT - 1 / aT));
    vaT = sqrt(mu * (2 / raT - 1 / aT));
    
    % DeltaV
    if type == 'pa'
        DeltaV1 = vpT - vp_i;
        DeltaV2 = va_f - vaT;
        om_f = om_i;
    elseif type == 'ap'
        DeltaV1 = vaT - va_i;
        DeltaV2 = vp_f - vpT;
        om_f = om_i;
    elseif type == 'pp'
        DeltaV1 = vpT - vp_i;
        DeltaV2 = vp_f - vaT;
        om_f = PositiveAngles(om_i + pi);
    elseif type == 'aa'
        DeltaV1 = vpT - va_i;
        DeltaV2 = va_f - vaT;

        if om_i + pi > 2*pi
            Deltat = 2 * pi * sqrt(aT^3 / mu);
        end
        
        om_f = PositiveAngles(om_i + pi);
    end
    
    % Tempo di trasferimento (metà periodo dell'orbita di trasferimento)
    Deltat = Deltat + pi * sqrt(aT^3 / mu);
end