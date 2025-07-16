function [dv1,dv2,dv3,thf,dt]=changeOrbitShape(ai,ei,omi,af,ef,omf,P,mu,ra_t)

% Change of Shape maneuver
%
% [dv1,dv2,dv3,thf,dt]=changeOrbitShape(ai,ei,omi,af,ef,omf,P,ra_t,mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai       [1x1]    initial semi-major axis                 [km]
% ei       [1x1]    initial eccentricity                    [-]
% omi      [1x1]    initial pericenter anomaly              [rad]
% af       [1x1]    final semi-major axis                   [km]
% ef       [1x1]    final eccentricity                      [-]
% omf      [1x1]    final pericenter anomaly                [rad]
% P - Maneuver point selection:
%       1: Perform maneuver at farthest point from current position
%       2: Perform maneuver at closest point to current position
% ra_t     [1x1]    manouver apocenter radius               [km]
% mu       [1x1]    gravitational parameter                 [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% dv1      [1x1]    first delta V                           [km/s]
% dv2      [1x1]    second delta V                          [rad]
% dv3      [1x1]    third delta V                           [rad]
% thf      [1x1]    final true anomaly after manouver       [rad]
% dt       [1x1]    manouver time                           [s]


% Rendo gli angoli positivi tra 0 e 2*pi
omf = PositiveAngles(omf);
omi = PositiveAngles(omi);

% Determina se il trasferimento è bitangente (8 argomenti) o biellittico (9)
isBitangent = (nargin == 8);

delta_omega = round(abs(omi - omf), 12);
isOpposite = (delta_omega == round(pi, 12));
isSame = (delta_omega == 0 || delta_omega == round(2*pi, 12));

% ---------------------------
% Caso BITANGENTE
% ---------------------------
if isBitangent
    % Determina tipo di trasferimento
    if P == 1 && isSame
        type = 'pa';
        thf = 0;
    elseif P == 2 && isSame
        type = 'ap';
        thf = 0;
    elseif P == 1 && isOpposite
        type = 'pp';
        thf = pi;
    elseif P == 2 && isOpposite
        type = 'aa';
        thf = pi;
    else
        error('Configurazione non supportata per trasferimento bitangente');
    end
    
    % Calcolo con funzione bitangentTransfer
    [dv1, dv2, dt, ~] = bitangentTransfer(ai, ei, af, ef, omi, type, mu);
    dv3 = 0;

% ---------------------------
% Caso BIELLITTICO
% ---------------------------
else
    % Calcola con funzione biellipticTransfer
    [dv1, dv2, dv3, dt1, dt2, ~] = biellipticTransfer(ai, ei, af, ef, ra_t, omi, omf, P, mu);
    dt = dt1 + dt2;
    
    % Imposta angolo finale thf
    if ((isSame && P == 1) || (isOpposite && P == 2))
        thf = 0;
    else
        thf = pi;
    end
end
end