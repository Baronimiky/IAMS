function [dv1,dv2,dv3,thf,dt,om_f_new]=changeOrbitShape(ai,ei,omi,af,ef,omf,P,ra_t)

    % thf è l'angolo finale sull'orbita di trasferimento
    
    mu=398600;

    % Rendi gli angoli positivi tra 0 e 2*pi
    omf = PositiveAngles(omf);
    omi = PositiveAngles(omi);

    % Determina se il trasferimento è bitangente (7 argomenti) o biellittico (8)
    isBitangent = (nargin == 7);

    % Calcola la differenza in omega
    delta_omega = round(abs(omi - omf), 12);
    
    angular_tol = 1e-10;
    isOpposite = abs(mod(delta_omega, 2*pi) - pi) < angular_tol;
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
        [dv1, dv2, dt, om_f_new] = bitangentTransfer(ai, ei, af, ef, omi, type, mu);
        dv3 = 0;

    % ---------------------------
    % Caso BIELLITTICO
    % ---------------------------
    else
        % Calcola con funzione biellipticTransfer
        [dv1, dv2, dv3, dt1, dt2, om_f_new] = biellipticTransfer(ai, ei, af, ef, ra_t, omi, mu);
        dt = dt1 + dt2;

        % Imposta angolo finale thf
        if (isSame && P == 1) || (isOpposite && P == 2)
            thf = 0;
        else
            thf = pi;
        end
    end
end