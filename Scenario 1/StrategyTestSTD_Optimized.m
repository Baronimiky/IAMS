function [optimal_time_strategies, optimal_deltav_strategies, minimum_time, minimum_deltav, ...
          all_times, all_deltav, strategy_parameters] = ...
          StrategyTestSTD_Optimized(ai, ei, ...
          initial_inclination, initial_raan, initial_periapsis_arg, initial_true_anomaly, ...
          final_semimajor_axis, final_eccentricity, final_inclination, final_raan, ...
          final_periapsis_arg, final_true_anomaly)

% OPTIMIZED_ORBITAL_TRANSFER_CALCULATOR - Calcola tutte le possibili strategie di trasferimento orbitale
%
% Questa funzione calcola 192 strategie standard di trasferimento orbitale, considerando
% tutte le possibili combinazioni di:
% - Ordine delle manovre: Piano-Assi-Forma (P-A-S), Piano-Forma-Assi (P-S-A), etc.
% - Tipi di trasferimento: bitangente e biellittico
% - Punti di manovra: pericentro/apocentro per ogni operazione
%
% -------------------------------------------------------------------------
% PARAMETRI DI INPUT:
% ai [1x1]    - Semiasse maggiore iniziale [km]
% ei [1x1]      - Eccentricità iniziale [-]
% initial_inclination [1x1]       - Inclinazione iniziale [rad]
% initial_raan [1x1]              - RAAN iniziale [rad]
% initial_periapsis_arg [1x1]     - Argomento del periapside iniziale [rad]
% initial_true_anomaly [1x1]      - Anomalia vera iniziale [rad]
% final_semimajor_axis [1x1]      - Semiasse maggiore finale [km]
% final_eccentricity [1x1]        - Eccentricità finale [-]
% final_inclination [1x1]         - Inclinazione finale [rad]
% final_raan [1x1]                - RAAN finale [rad]
% final_periapsis_arg [1x1]       - Argomento del periapside finale [rad]
% final_true_anomaly [1x1]        - Anomalia vera finale [rad]
%
% -------------------------------------------------------------------------
% PARAMETRI DI OUTPUT:
% optimal_time_strategies [Nx5]    - Strategie ottimali per tempo minimo
%                                   [sequence_code, transfer_type, plane_point, periapsis_point, deltav_cost]
% optimal_deltav_strategies [Nx5]  - Strategie ottimali per delta-v minimo
%                                   [sequence_code, transfer_type, plane_point, periapsis_point, time_cost]
% minimum_time [1x1]               - Tempo minimo di trasferimento [s]
% minimum_deltav [1x1]             - Delta-v minimo di trasferimento [km/s]
% all_times [192x1]                - Vettore di tutti i tempi calcolati [s]
% all_deltav [192x1]               - Vettore di tutti i delta-v calcolati [km/s]
% strategy_parameters [192x17]     - Parametri dettagliati di tutte le strategie
%                                   [bielliptical_radius, deltav_plane, deltav_periapsis, deltav1, deltav2, deltav3,
%                                    time_plane, time_periapsis, time_transit, time_shape, time_final,
%                                    final_periapsis_arg_plane, final_periapsis_arg_periapsis, maneuver_true_anomaly_plane,
%                                    old_final_true_anomaly, new_final_true_anomaly, shape_true_anomaly]
%
% -------------------------------------------------------------------------
% LEGENDA CODICI NUMERICI:
% sequence_code (1-6):
%   1: Piano-Assi-Forma    2: Piano-Forma-Assi    3: Assi-Piano-Forma
%   4: Forma-Piano-Assi    5: Assi-Forma-Piano    6: Forma-Assi-Piano
%
% transfer_type (1-8):
%   1: bitangente Pericentro-Apocentro     2: bitangente Pericentro-Pericentro
%   3: bitangente Apocentro-Pericentro     4: bitangente Apocentro-Apocentro
%   5: biellittico Pericentro-Pericentro   6: biellittico Pericentro-Apocentro
%   7: biellittico Apocentro-Apocentro     8: biellittico Apocentro-Pericentro
%
% plane_point, periapsis_point (1-2): Punto di manovra (1 o 2)

% Parametro gravitazionale standard
mu = 398600;

% Numero totale di strategie da calcolare
total_strategies = 192;

% Pre-allocazione delle matricri di output per efficienza
strategy_codes = zeros(total_strategies, 4);
all_times = zeros(total_strategies, 1);
all_deltav = zeros(total_strategies, 1);
strategy_parameters = zeros(total_strategies, 17);

% Pre-calcolo dei parametri per trasferimenti biellittici
final_radius_periapsis = final_semimajor_axis * (1 - final_eccentricity^2) / (1 - final_eccentricity);
bielliptical_radius_range = 0:0.1:5*final_radius_periapsis;
num_bielliptical_radii = length(bielliptical_radius_range);

% Cache per evitare calcoli ripetuti
plane_change_cache = containers.Map();
shape_change_cache = containers.Map();
time_flight_cache = containers.Map();
bielliptical_cache = containers.Map();

% Pre-calcolo dei parametri dipendenti dal tipo di trasferimento
phase_shift_values = [0, pi, 0, pi, 0, pi, 0, pi];  % per transfer_type = 1:8
transfer_point_values = [1, 1, 2, 2, 1, 1, 2, 2];  % per transfer_type = 1:8
shape_true_anomaly_values = [2*pi, 2*pi, pi, pi, 2*pi, 2*pi, pi, pi];  % per transfer_type = 1:8

strategy_counter = 1;

% Loop principale attraverso tutte le combinazioni di strategie
for sequence_code = 1:6
    for transfer_type = 1:8
        phase_shift = phase_shift_values(transfer_type);
        transfer_point = transfer_point_values(transfer_type);
        shape_true_anomaly = shape_true_anomaly_values(transfer_type);
        
        % Pre-calcolo per trasferimenti biellittici
        if transfer_type > 4
            bielliptical_key = sprintf('type_%d', transfer_type);
            if ~isKey(bielliptical_cache, bielliptical_key)
                [bielliptical_deltav_vec, bielliptical_time_vec, bielliptical_omf_vec] = ...
                    compute_bielliptical_options(bielliptical_radius_range, ai, ...
                    ei, final_semimajor_axis, final_eccentricity, transfer_point);
                bielliptical_cache(bielliptical_key) = struct('deltav', bielliptical_deltav_vec, ...
                    'time', bielliptical_time_vec, 'omf', bielliptical_omf_vec);
            end
            cached_bielliptical = bielliptical_cache(bielliptical_key);
        end
        
        for plane_maneuver_point = 1:2
            for periapsis_maneuver_point = 1:2
                bielliptical_radius = 0;
                
                % Calcolo specifico per ogni sequenza di manovre
                switch sequence_code
                    case 1 % Piano-Assi-Forma
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(ai, ei, initial_inclination, ...
                            initial_raan, initial_periapsis_arg, final_inclination, final_raan, ...
                            initial_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(ai, ei, plane_final_periapsis_arg, ...
                            final_periapsis_arg + phase_shift, plane_true_anomaly, periapsis_maneuver_point,mu);
                        
                        transit_time = cached_time_of_flight(ai, ei, ...
                            new_final_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, periapsis_final_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, final_periapsis_arg, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            shape_final_true_anomaly, final_true_anomaly, time_flight_cache, mu);
                        total_shape_time = shape_time + transit_time;
                        
                    case 2 % Piano-Forma-Assi
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(ai, ei, initial_inclination, ...
                            initial_raan, initial_periapsis_arg, final_inclination, final_raan, ...
                            initial_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        
                        transit_time = cached_time_of_flight(ai, ei, ...
                            plane_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, plane_final_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, plane_final_periapsis_arg + phase_shift, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        total_shape_time = shape_time + transit_time;
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(final_semimajor_axis, final_eccentricity, shape_final_periapsis_arg, ...
                            final_periapsis_arg, plane_true_anomaly, periapsis_maneuver_point,mu);
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            new_final_true_anomaly, final_true_anomaly, time_flight_cache , mu);
                        
                    case 3 % Assi-Piano-Forma
                        temp_periapsis_arg = InverseCP(initial_inclination, initial_raan, final_periapsis_arg - phase_shift, final_inclination, final_raan);
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(ai, ei, initial_periapsis_arg, ...
                            temp_periapsis_arg, initial_true_anomaly, periapsis_maneuver_point,mu);
                        
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(ai, ei, initial_inclination, ...
                            initial_raan, periapsis_final_periapsis_arg, final_inclination, final_raan, ...
                            new_final_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        
                        transit_time = cached_time_of_flight(ai, ei, ...
                            plane_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, plane_final_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, final_periapsis_arg, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        total_shape_time = shape_time + transit_time;
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            shape_final_true_anomaly, final_true_anomaly, time_flight_cache, mu);
                        
                    case 4 % Forma-Piano-Assi
                        transit_time = cached_time_of_flight(ai, ei, ...
                            initial_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, initial_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, initial_periapsis_arg + phase_shift, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        total_shape_time = shape_time + transit_time;
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(final_semimajor_axis, final_eccentricity, initial_inclination, ...
                            initial_raan, shape_final_periapsis_arg, final_inclination, final_raan, ...
                            shape_final_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(final_semimajor_axis, final_eccentricity, plane_final_periapsis_arg, ...
                            final_periapsis_arg, plane_true_anomaly, periapsis_maneuver_point,mu);
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            new_final_true_anomaly, final_true_anomaly, time_flight_cache, mu);
                        
                    case 5 % Assi-Forma-Piano
                        temp_periapsis_arg = InverseCP(initial_inclination, initial_raan, final_periapsis_arg, final_inclination, final_raan);
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(ai, ei, initial_periapsis_arg, ...
                            temp_periapsis_arg + phase_shift, initial_true_anomaly, periapsis_maneuver_point,mu);
                        
                        transit_time = cached_time_of_flight(ai, ei, ...
                            new_final_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, periapsis_final_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, temp_periapsis_arg, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        total_shape_time = shape_time + transit_time;
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(final_semimajor_axis, final_eccentricity, initial_inclination, ...
                            initial_raan, temp_periapsis_arg, final_inclination, final_raan, ...
                            shape_final_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            plane_true_anomaly, final_true_anomaly, time_flight_cache, mu);
                        
                    case 6 % Forma-Assi-Piano
                        transit_time = cached_time_of_flight(ai, ei, ...
                            initial_true_anomaly, shape_true_anomaly, time_flight_cache, mu);
                        
                        if transfer_type > 4
                            [shape_deltav, shape_time, bielliptical_radius, shape_final_periapsis_arg] = ...
                                get_optimal_bielliptical(cached_bielliptical, bielliptical_radius_range);
                            shape_deltav1 = 0; shape_deltav2 = 0; shape_deltav3 = 0;
                            shape_final_true_anomaly = shape_true_anomaly;
                        else
                            [shape_deltav1, shape_deltav2, shape_deltav3, shape_final_true_anomaly, shape_time, shape_final_periapsis_arg] = ...
                                cached_shape_change(ai, ei, initial_periapsis_arg, ...
                                final_semimajor_axis, final_eccentricity, initial_periapsis_arg + phase_shift, transfer_point, shape_change_cache);
                            shape_deltav = abs(shape_deltav1) + abs(shape_deltav2) + abs(shape_deltav3);
                        end
                        
                        total_shape_time = shape_time + transit_time;
                        temp_periapsis_arg = InverseCP(initial_inclination, initial_raan, final_periapsis_arg, final_inclination, final_raan);
                        [periapsis_deltav, periapsis_final_periapsis_arg, new_final_true_anomaly, old_final_true_anomaly, periapsis_time] = ...
                            changePericenterArg(final_semimajor_axis, final_eccentricity, shape_final_periapsis_arg, ...
                            temp_periapsis_arg, shape_final_true_anomaly, periapsis_maneuver_point,mu);
                        
                        [plane_deltav, plane_final_periapsis_arg, plane_true_anomaly, plane_time] = ...
                            cached_plane_change(final_semimajor_axis, final_eccentricity, initial_inclination, ...
                            initial_raan, temp_periapsis_arg, final_inclination, final_raan, ...
                            new_final_true_anomaly, plane_maneuver_point, mu, plane_change_cache);
                        final_time = cached_time_of_flight(final_semimajor_axis, final_eccentricity, ...
                            plane_true_anomaly, final_true_anomaly, time_flight_cache, mu);
                end
                
                % Memorizzazione dei risultati
                strategy_parameters(strategy_counter,:) = [bielliptical_radius, plane_deltav, periapsis_deltav, ...
                    shape_deltav1, shape_deltav2, shape_deltav3, plane_time, periapsis_time, transit_time, ...
                    shape_time, final_time, plane_final_periapsis_arg, periapsis_final_periapsis_arg, ...
                    plane_true_anomaly, old_final_true_anomaly, new_final_true_anomaly, shape_final_true_anomaly];
                
                all_deltav(strategy_counter) = abs(shape_deltav) + abs(plane_deltav) + abs(periapsis_deltav);
                all_times(strategy_counter) = plane_time + total_shape_time + periapsis_time + final_time;
                strategy_codes(strategy_counter,:) = [sequence_code, transfer_type, plane_maneuver_point, periapsis_maneuver_point];
                strategy_counter = strategy_counter + 1;
            end
        end
    end
end

% Ricerca dei valori ottimali usando operazioni vettoriali
[minimum_deltav, ~] = min(all_deltav);
[minimum_time, ~] = min(all_times);

% Identificazione delle strategie ottimali
optimal_deltav_indices = find(all_deltav == minimum_deltav);
optimal_time_indices = find(all_times == minimum_time);

% Costruzione delle matrici di output
optimal_deltav_strategies = [strategy_codes(optimal_deltav_indices, :), all_times(optimal_deltav_indices)];
optimal_time_strategies = [strategy_codes(optimal_time_indices, :), all_deltav(optimal_time_indices)];

end

% -------------------------------------------------------------------------
% FUNZIONI AUSILIARIE PER IL CACHING E L'OTTIMIZZAZIONE
% -------------------------------------------------------------------------

function [bielliptical_deltav_vec, bielliptical_time_vec, bielliptical_omf_vec] = compute_bielliptical_options(radius_range, ai, ei, af, ef, transfer_point)
    % Calcola tutte le opzioni per i trasferimenti biellittici
    num_radii = length(radius_range);
    bielliptical_deltav_vec = zeros(1, num_radii);
    bielliptical_time_vec = zeros(1, num_radii);
    bielliptical_omf_vec = zeros(1, num_radii);
    
    for i = 1:num_radii
        [deltav1, deltav2, deltav3, ~, time_temp] = changeOrbitShape(ai, ei, 0, af, ef, 0, transfer_point, radius_range(i));
        bielliptical_deltav_vec(i) = abs(deltav1) + abs(deltav2) + abs(deltav3);
        bielliptical_time_vec(i) = time_temp;
        % Per i trasferimenti biellittici, l'omf finale è tipicamente uguale a quello iniziale
        bielliptical_omf_vec(i) = 0;
    end
end

function [optimal_deltav, optimal_time, optimal_radius, optimal_omf] = get_optimal_bielliptical(cached_data, radius_range)
    % Estrae la soluzione ottimale dai dati biellittici pre-calcolati
    [optimal_deltav, min_index] = min(cached_data.deltav);
    optimal_time = cached_data.time(min_index);
    optimal_radius = radius_range(min_index);
    optimal_omf = cached_data.omf(min_index);
end

function [deltav, final_periapsis_arg, true_anomaly, time] = cached_plane_change(a, e, initial_inclination, initial_raan, initial_periapsis_arg, final_inclination, final_raan, initial_true_anomaly, maneuver_point, mu, cache)
    % Funzione con cache per il cambio di piano orbitale
    cache_key = sprintf('%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%d', ...
        a, e, initial_inclination, initial_raan, initial_periapsis_arg, final_inclination, final_raan, initial_true_anomaly, maneuver_point);
    
    if isKey(cache, cache_key)
        result = cache(cache_key);
        deltav = result.deltav;
        final_periapsis_arg = result.final_periapsis_arg;
        true_anomaly = result.true_anomaly;
        time = result.time;
    else
        [deltav, final_periapsis_arg, true_anomaly] = changeOrbitalPlane(a, e, initial_inclination, initial_raan, initial_periapsis_arg, final_inclination, final_raan, initial_true_anomaly, mu);
        time = 0; % Assumendo che il cambio di piano sia istantaneo
        cache(cache_key) = struct('deltav', deltav, 'final_periapsis_arg', final_periapsis_arg, 'true_anomaly', true_anomaly, 'time', time);
    end
end

function [deltav1, deltav2, deltav3, final_true_anomaly, time, final_periapsis_arg] = ...
    cached_shape_change(ai, ei, initial_periapsis_arg, af, ef, final_periapsis_arg, transfer_point, cache)

    cache_key = sprintf('%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%d', ...
        ai, ei, initial_periapsis_arg, af, ef, final_periapsis_arg, transfer_point);

    if isKey(cache, cache_key)
        result = cache(cache_key);
        deltav1 = result.deltav1;
        deltav2 = result.deltav2;
        deltav3 = result.deltav3;
        final_true_anomaly = result.final_true_anomaly;
        time = result.time;
        final_periapsis_arg = result.final_periapsis_arg;
    else
        [deltav1, deltav2, deltav3, final_true_anomaly, time, final_periapsis_arg] = ...
            changeOrbitShape(ai, ei, initial_periapsis_arg, af, ef, final_periapsis_arg, transfer_point);
        
        cache(cache_key) = struct( ...
            'deltav1', deltav1, ...
            'deltav2', deltav2, ...
            'deltav3', deltav3, ...
            'final_true_anomaly', final_true_anomaly, ...
            'time', time, ...
            'final_periapsis_arg', final_periapsis_arg);
    end
end


function time = cached_time_of_flight(a, e, initial_true_anomaly, final_true_anomaly, cache, mu)
    % Funzione con cache per il calcolo del tempo di volo
    cache_key = sprintf('%.6f_%.6f_%.6f_%.6f', a, e, initial_true_anomaly, final_true_anomaly);
    
    if isKey(cache, cache_key)
        time = cache(cache_key);
    else
        time = TOF(a, e, initial_true_anomaly, final_true_anomaly,mu);
        cache(cache_key) = time;
    end
end