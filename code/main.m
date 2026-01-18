function main()
    % PROJECT_ANALYSIS
    % File unico per l'analisi del sistema dei serbatoi in cascata.
    % Include: Parametri, Equilibrio, Linearizzazione, FdT e Bode.

    clc;
    close all; % Chiude eventuali figure aperte
    fprintf('--------------------------------------------------\n');
    fprintf('PROGETTO CONTROLLI AUTOMATICI - GRUPPO 35\n');
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
    u_eq = (k1 * sqrt(x1_eq_target)) / k4;
    
    % Vettore di stato all'equilibrio
    x_eq = [x1_eq_target; x2_eq_target];
    
    %% 3. VERIFICA DELL'EQUILIBRIO
    dot_x1 = -k1 * sqrt(x_eq(1)) + k4 * u_eq;
    dot_x2 = k2 * sqrt(x_eq(1)) - k3 * sqrt(x_eq(2));

    fprintf('>>> 1. PUNTO DI EQUILIBRIO <<<\n');
    fprintf('x_eq = [%.4f; %.4f]\n', x_eq(1), x_eq(2));
    fprintf('u_eq = %.5f V\n', u_eq);
    
    if abs(dot_x1) < 1e-6 && abs(dot_x2) < 1e-3
        fprintf('[OK] Equilibrio verificato (derivate ~ 0).\n\n');
    else
        fprintf('[WARNING] Residui alti nelle derivate.\n\n');
    end

    %% 4. LINEARIZZAZIONE (Calcolo matrici A, B, C, D)
    % Calcoliamo i termini delle Jacobiane valutate in x_eq
    % Termini ricorrenti
    sqrt_x1 = sqrt(x_eq(1));
    sqrt_x2 = sqrt(x_eq(2));

    % Matrice A: df/dx
    % A11 = -k1 / (2*sqrt(x1))
    % A21 = k2 / (2*sqrt(x1))
    % A22 = -k3 / (2*sqrt(x2))
    A = [ -k1/(2*sqrt_x1),       0;
           k2/(2*sqrt_x1),  -k3/(2*sqrt_x2) ];

    % Matrice B: df/du
    B = [ k4; 
          0 ];

    % Matrici C e D (y = x2)
    C = [0, 1];
    D = 0;

    fprintf('>>> 2. SISTEMA LINEARIZZATO (Matrici) <<<\n');
    disp('Matrice A:'); disp(A);
    disp('Matrice B:'); disp(B);
    disp('Matrice C:'); disp(C);
    
    % Analisi autovalori (poli)
    poli = eig(A);
    fprintf('Poli del sistema (autovalori di A): \n');
    disp(poli);

    %% 5. FUNZIONE DI TRASFERIMENTO (Punto 2)
    % Creiamo l'oggetto state-space (ss) e lo convertiamo in transfer function (tf)
    sys_ss = ss(A, B, C, D);
    sys_tf = tf(sys_ss);

    fprintf('>>> 3. FUNZIONE DI TRASFERIMENTO G(s) <<<\n');
    % Stampa a video la G(s)
    sys_tf 
    
    % Verifica guadagno statico
    mu = dcgain(sys_tf);
    fprintf('Guadagno statico mu: %.4f (%.2f dB)\n\n', mu, 20*log10(mu));

    %% 6. DIAGRAMMA DI BODE E SALVATAGGIO
    fprintf('>>> 4. GENERAZIONE GRAFICI <<<\n');
    
    % Setup della figura
    fig_bode = figure('Name', 'Diagramma di Bode', 'Color', 'w');
    opts = bodeoptions;
    opts.FreqUnits = 'rad/s';
    opts.Grid = 'on';
    opts.Title.String = 'Diagramma di Bode - G(s)';
    opts.PhaseMatching = 'on'; % Per evitare salti di fase strani
    
    % Plot
    bode(sys_tf, opts);
    
    % Personalizzazione linee (più spesse per il report)
    lines = findall(gcf,'type','line');
    set(lines, 'LineWidth', 1.5);

    % Gestione salvataggio nella cartella 'figs'
    folder_name = 'figs';
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
        fprintf('Cartella "%s" creata.\n', folder_name);
    end
    
    filename = fullfile(folder_name, 'bode.jpg');
    
    % Salvataggio ad alta risoluzione (300 DPI)
    print(fig_bode, filename, '-djpeg', '-r300');
    fprintf('Immagine salvata correttamente in: %s\n', filename);
    
    fprintf('--------------------------------------------------\n');
    fprintf('ANALISI COMPLETATA.\n');
end