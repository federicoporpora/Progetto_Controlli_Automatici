function main()
    % PROJECT_ANALYSIS
    % File unico per l'analisi del sistema dei serbatoi in cascata.
    % Esegue i calcoli dei parametri e dell'equilibrio senza salvare
    % variabili nel workspace globale.

    clc;
    fprintf('--------------------------------------------------\n');
    fprintf('PROGETTO CONTROLLI AUTOMATICI - ANALISI SISTEMA\n');
    fprintf('--------------------------------------------------\n\n');

    %% 1. DEFINIZIONE PARAMETRI
    % Dati tratti dalla Tabella 1 della traccia
    k1 = 0.10;      % Parametro scarico serbatoio 1
    k2 = 1.50;      % Parametro ingresso serbatoio 2
    k3 = 0.12;      % Parametro scarico serbatoio 2
    k4 = 2.50;      % Guadagno pompa
    
    % Valori di equilibrio desiderati (da Tabella 1)
    x1_eq_target = 0.0235; % Livello a1 all'equilibrio [m]
    x2_eq_target = 3.67;   % Livello a2 all'equilibrio [m]

    %% 2. CALCOLO DELL'EQUILIBRIO (Punto 1)
    % Risolviamo il sistema f(x,u) = 0 per trovare u_eq.
    % Eq 1: -k1*sqrt(x1) + k4*u = 0  => u = (k1*sqrt(x1))/k4
    
    u_eq = (k1 * sqrt(x1_eq_target)) / k4;
    
    % Vettore di stato all'equilibrio
    x_eq = [x1_eq_target; x2_eq_target];
    
    %% 3. VERIFICA DELL'EQUILIBRIO
    % Sostituiamo i valori calcolati nelle equazioni differenziali originali
    % per assicurarci che le derivate siano nulle (o quasi nulle).
    
    % dx1 = -k1*sqrt(x1) + k4*u
    dot_x1 = -k1 * sqrt(x_eq(1)) + k4 * u_eq;
    
    % dx2 = k2*sqrt(x1) - k3*sqrt(x2)
    dot_x2 = k2 * sqrt(x_eq(1)) - k3 * sqrt(x_eq(2));

    %% 4. STAMPA DEI RISULTATI
    fprintf('>>> PUNTO DI EQUILIBRIO CALCOLATO <<<\n');
    fprintf('Stato x1 (Livello T1): %10.4f [m]\n', x_eq(1));
    fprintf('Stato x2 (Livello T2): %10.4f [m]\n', x_eq(2));
    fprintf('Ingresso u (Tensione): %10.5f [V]\n', u_eq);
    fprintf('\n');
    
    fprintf('>>> VERIFICA DERIVATE (Devono essere ~0) <<<\n');
    fprintf('dx1/dt: %10.4e\n', dot_x1);
    fprintf('dx2/dt: %10.4e\n', dot_x2);
    
    if abs(dot_x1) < 1e-6 && abs(dot_x2) < 1e-4
        fprintf('\n[OK] Il punto calcolato è un equilibrio valido.\n');
    else
        fprintf('\n[WARNING] I residui sono alti. Controllare i parametri.\n');
        % Nota: dx2 potrebbe non essere esattamente 0 se i dati 
        % forniti in tabella (3.67) sono arrotondati.
    end

end
