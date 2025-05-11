function deltat = TOF(a,e,th1,th2,mu)

%TOF Time Of Flight

%Input arg:
% a   [1x1]   semi-major axis         [km]
% e   [1x1]   eccentricity            [-]
% th1 [1x1]   Initial true anomaly    [rad]
% th2 [1x1]   Final true anomaly      [rad]
% mu  [1x1]   gravitational paramete  [km^3/s^2]
%----------------------------------------------
%Output arg:
% deltat   [1x1]   time of flight       [s]
%
%-----------------------------------------------

%anomalie eccentriche

T = 2*pi*sqrt((abs(a).^3)/mu); % periodo utile solo per orbite ellittiche

%
% Nota:
%   Se theta2 risulta minore di theta1, si assume che theta2 appartenga al giro
%   successivo e si aggiunge 2*pi. La conversione in anomalia eccentrica
%   utilizza atan2 per evitare problemi numerici vicino a 0 e 2*pi.

    % Se theta2 è minore di theta1, si assume che theta2 appartenga al giro successivo.
    if th2 < th1
        th2 = th2 + 2*pi;
    end

    if e >= 0 && e < 1
        % orbita ellittica: conversione stabile
        E1 = 2 * atan2( sqrt(1 - e^2) * sin(th1), 1 + e * cos(th1) );
        E2 = 2 * atan2( sqrt(1 - e^2) * sin(th2), 1 + e * cos(th2) );

        if E2 < E1
            E2 = E2 + 2*pi;
        end

        M1 = E1 - e*sin(E1);
        M2 = E2 - e*sin(E2);

        deltat = sqrt(a^3 / mu) * (M2 - M1)

    elseif e > 1
        % orbita iperbolica
        arg1 = tan(th1/2) * sqrt((e-1)/(e+1));
        arg2 = tan(th2/2) * sqrt((e-1)/(e+1));
        
        % Protezione: limita argomenti entro dominio valido [-0.9999, 0.9999]
        arg1 = min(max(arg1, -0.9999), 0.9999);
        arg2 = min(max(arg2, -0.9999), 0.9999);
        
        F1 = 2 * atanh(arg1);
        F2 = 2 * atanh(arg2);

        M1 = e*sinh(F1) - F1;
        M2 = e*sinh(F2) - F2;

        a_abs = abs(a); % a è negativo per le orbite iperboliche
        deltat = abs(sqrt(a_abs^3 / mu) * (M2 - M1))

    else
        error('Orbita parabolica (e=1) non supportata in questa funzione.');
    end
end