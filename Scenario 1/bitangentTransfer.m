function [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(ai, ei, af, ef, type, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(ai, ei, af, ef, type, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai        [1x1]  initial semi-major axis          [km]
% ei        [1x1]  initial eccentricity             [-]
% af        [1x1]  final semi-major axis            [km]
% ef        [1x1]  final eccentricity               [-]
% type      [char]  maneuver type   ('pa', 'ap', 'pp', 'aa')
% mu        [1x1]  gravitational parameter          [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1   [1x1]  1st maneuver impulse            [km/s]
% DeltaV2   [1x1]  2nd maneuver impulse            [km/s]
% Deltat    [1x1]  maneuver time                   [s]
%

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
            raT = ra_i;
        case 'aa'  % Apocentro a Apocentro
            rpT = ra_i;
            raT = ra_f;
        otherwise
            error('Invalid maneuver type. Choose pa, ap, pp, or aa.');
    end
    
    % Compute transfer orbit semi-major axis
    aT = (rpT + raT) / 2;
    
    % Compute velocity values
    vp_i = sqrt(mu * (2 / rp_i - 1 / ai));
    va_i = sqrt(mu * (2 / ra_i - 1 / ai));
    vp_f = sqrt(mu * (2 / rp_f - 1 / af));
    va_f = sqrt(mu * (2 / ra_f - 1 / af));
    vpT = sqrt(mu * (2 / rpT - 1 / aT));
    vaT = sqrt(mu * (2 / raT - 1 / aT));
    
    % Compute DeltaV values
    if type == "pa"
        DeltaV1 = vpT - vp_i;
        DeltaV2 = va_f - vaT;
    elseif type == "ap"
        DeltaV1 = vaT - va_i;
        DeltaV2 = vp_f - vpT;
    elseif type == "pp"
        DeltaV1 = vpT - vp_i;
        DeltaV2 = vp_f - vpT;
    elseif type == "aa"
        DeltaV1 = vaT - va_i;
        DeltaV2 = va_f - vaT;
    end
    
    % Compute transfer time (half period of transfer orbit)
    Deltat = pi * sqrt(aT^3 / mu);
end
