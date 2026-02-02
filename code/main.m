close all; clear all; clc;

% Setup Cartella Output
output_dir = fullfile('..', 'report', 'figs');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% PARAMETRI DEL SISTEMA (Serbatoi) %%
k1 = 0.10;
k2 = 1.50;
k3 = 0.12;
k4 = 2.50;

% Valori di equilibrio (da traccia)
a1_eq = 0.0235;
a2_eq = 3.67;

%% PUNTO 1 %%
% Calcolo dell'equilibrio e linearizzazione
x1e = a1_eq;
x2e = a2_eq;
xe = [x1e; x2e];

% Ingresso di equilibrio (Tensione stazionaria della pompa)
ue = (k1 * sqrt(x1e)) / k4;
fprintf('Tensione di equilibrio ue = %.4f V\n', ue);

% Matrici del sistema linearizzato (Jacobiane)
A = [ -k1/(2*sqrt(x1e)),      0;
       k2/(2*sqrt(x1e)), -k3/(2*sqrt(x2e)) ];
B = [ k4; 0 ];
C = [ 0, 1 ];
D = 0;

fprintf('PUNTO 1: Equilibrio calcolato.\n');

%% PUNTO 2 %%
% Calcolo Funzione di Trasferimento G(s)
sys_ss = ss(A, B, C, D);
s = tf('s');

% Calcolo G(s) analitico
poli = eig(A);
Num_val = B(1) * A(2,1); 
G = Num_val / ((s - poli(1)) * (s - poli(2)));

% Stampa Grafico Bode G(s)
figure(1); clf; set(gcf, 'Color', 'w');
margin(G);
title('Diagramma di Bode di G(s)');
grid on;
saveas(gcf, fullfile(output_dir, 'bode_G_plant.jpg'));

%% PUNTO 3 %%
% SINTESI DEL REGOLATORE

% Specifiche di Progetto
e_star = 0.05;      
W_star = 3;         
D_star = 2.5;       
Mf_target = 55;     
wc_des = 160;       % Pulsazione critica target (rad/s)

% -- Regolatore Statico --
mu_G = dcgain(G);
Kr = 250; 
R_s = Kr;
G_e = G * R_s; 

% -- Regolatore Dinamico (Rete Anticipatrice) --
[mag_Ge, phase_Ge] = bode(G_e, wc_des);
val_Ge = mag_Ge * exp(1j*deg2rad(phase_Ge));

M_star = 1 / abs(val_Ge);
phi_star_deg = Mf_target - (180 + rad2deg(angle(val_Ge)));

while phi_star_deg > 85
    phi_star_deg = phi_star_deg - 5; 
end

phi_rad = deg2rad(phi_star_deg);
tau = (M_star - cos(phi_rad)) / (wc_des * sin(phi_rad));
alpha_tau = (cos(phi_rad) - 1/M_star) / (wc_des * sin(phi_rad));
alpha = alpha_tau / tau;

% Costruzione Regolatore
R_d = (1 + tau*s) / (1 + alpha*tau*s);
R_tot = R_s * R_d; % Regolatore Totale
L = G * R_tot;     % Funzione d'anello

% -- CALCOLO DATI ESATTI PER IL GRAFICO BODE --
[Gm, Pm, Wcg, Wcp] = margin(L); 
phase_at_Wcp = -180 + Pm; 
limit_phase = -150; 

w_vec = logspace(-2, 7, 500);
[mag_L, phase_L] = bode(L, w_vec);
mag_L_db = 20*log10(squeeze(mag_L));
phase_L_deg = squeeze(phase_L);

% -- GRAFICO LOOP SHAPING --
figure(2); clf; set(gcf, 'Color', 'w');

% Subplot 1: MODULO (dB)
subplot(2,1,1);
semilogx(w_vec, mag_L_db, 'b', 'LineWidth', 1.5); grid on; hold on;
ylabel('Magnitude (dB)'); title('Loop Shaping L(s) - Modulo');
xlim([1e-2 1e7]); ylim([-150 100]);

patch([1e-4 2 2 1e-4], [-200 -200 40 40], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r');
text(0.05, 20, 'Disturbi (>40dB)', 'Color', 'r', 'FontWeight', 'bold');
yline(40, 'r--', 'HandleVisibility', 'off');

patch([1e5 1e7 1e7 1e5], [-63 -63 150 150], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r');
text(2e5, -40, 'Rumore (<-63dB)', 'Color', 'r', 'FontWeight', 'bold');
yline(-63, 'r--', 'HandleVisibility', 'off');

yline(0, 'k-');
plot(Wcp, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
text(Wcp*1.5, 10, sprintf('\\omega_c = %.1f rad/s', Wcp), 'Color', 'k');
xline(Wcp, 'k--');

% Subplot 2: FASE (Gradi)
subplot(2,1,2);
semilogx(w_vec, phase_L_deg, 'b', 'LineWidth', 1.5); grid on; hold on;
ylabel('Phase (deg)'); xlabel('Frequency (rad/s)');
title('Loop Shaping L(s) - Fase');
xlim([1e-2 1e7]); ylim([-200 -50]);
yline(-180, 'k-', 'Limit -180');
xline(Wcp, 'k--');

plot(Wcp, limit_phase, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
text(Wcp*1.3, limit_phase, ' Limite Minimo (-150^\circ)', 'Color', 'r', 'FontSize', 9);
plot(Wcp, phase_at_Wcp, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
text(Wcp*1.3, phase_at_Wcp, sprintf(' Fase Attuale: %.1f^\\circ', phase_at_Wcp), 'Color', 'g', 'FontWeight', 'bold');
quiver(Wcp, limit_phase-15, 0, 10, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');

saveas(gcf, fullfile(output_dir, 'loop_shaping.jpg'));

% -- Grafico Bode Comparativo --
figure(3); clf; set(gcf, 'Color', 'w');
bodemag(G, 'b--'); hold on;
bodemag(L, 'r-');
grid on;
legend('G(s) Plant', 'L(s) Open Loop');
title('Confronto G(s) vs L(s)');
saveas(gcf, fullfile(output_dir, 'bode_comparativo.jpg'));

%% PUNTO 4 %%
% TEST E SIMULAZIONI
S = 1 / (1 + L); 
F = L / (1 + L); 
Q_u = R_tot / (1 + L); % Funzione di Sensitività del Controllo (Ref -> Input u)

% -- Risposta Gradino W=3 (Verifica Completa Dinamica + Statica) --
figure(4); clf; set(gcf, 'Color', 'w');
step(W_star * F, 0.2);
grid on; hold on;

% Vincoli grafici
yline(W_star, 'k-', 'Riferimento');
yline(W_star * 1.05, 'g--', 'Settling 5%');
yline(W_star * 0.95, 'g--', 'HandleVisibility', 'off');
yline(W_star * 1.20, 'm--', 'Max Overshoot 20%');
yline(W_star - e_star, 'r-.', 'Errore Regime');
yline(W_star + e_star, 'r-.', 'HandleVisibility', 'off');
xline(0.050, 'k-', 'Max Time (0.05s)');

title(['Risposta al gradino W=' num2str(W_star) ' - Uscita y(t)']);
legend('Risposta', 'Riferimento', 'Settling 5%', 'Max Overshoot', 'Err. Regime', 'Max Time', 'Location', 'SouthEast');
saveas(gcf, fullfile(output_dir, 'step_response_output.jpg'));


% -- NUOVO: Analisi Sforzo di Controllo (Variabile Manipolata) --
figure(8); clf; set(gcf, 'Color', 'w');
% Risposta al gradino del segnale di controllo (delta u)
[u_step, t_u] = step(W_star * Q_u, 0.2); 
plot(t_u, u_step, 'b', 'LineWidth', 1.5);
grid on; hold on;
yline(0, 'k-');
title(['Sforzo di Controllo (Variazione \deltau) per W=' num2str(W_star)]);
ylabel('Tensione \deltau (Volt)'); xlabel('Tempo (s)');
max_u = max(abs(u_step));
text(0.1, max_u*0.8, sprintf('Max \\DeltaV \\approx %.1f V', max_u), 'BackgroundColor', 'w');
saveas(gcf, fullfile(output_dir, 'step_response_control_effort.jpg'));
fprintf('Picco sforzo di controllo (delta): %.2f V\n', max_u);


% -- Funzioni di Sensitività --
figure(5); clf; set(gcf, 'Color', 'w');
bodemag(S, 'b'); hold on;
bodemag(F, 'r--');
grid on;
legend('|S|', '|F|');
title('Funzioni di Sensitività S e F');
saveas(gcf, fullfile(output_dir, 'sensitivita.jpg'));

% -- SIMULAZIONE DISTURBI E RUMORE --

% 1. Disturbo sull'uscita d(t) -- MODIFICATO: Tempo esteso per regime
figure(6); clf; set(gcf, 'Color', 'w');

td = 0:0.1:150; % Esteso a 150s per garantire regime
d_val = zeros(size(td));
for k = 1:4
    d_val = d_val + 1.0 * sin(0.1 * k * td);
end

sim_d = lsim(S, d_val, td);

subplot(2,1,1);
plot(td, d_val, 'r--', 'LineWidth', 1);
ylabel('Ampiezza'); title('Ingresso: Disturbo d(t)');
grid on; xlim([0 max(td)]);

subplot(2,1,2);
plot(td, sim_d, 'b', 'LineWidth', 1.5);
ylabel('Ampiezza'); xlabel('Tempo (s)');
title('Uscita: Residuo y_d(t) (Attenuato)');
grid on; xlim([0 max(td)]);

% Calcolo attenuazione su seconda metà del segnale (REGIME)
idx_regime = round(length(td)/2):length(td);
max_in_d = max(abs(d_val(idx_regime)));
max_out_d = max(abs(sim_d(idx_regime)));
att_d_db = 20*log10(max_in_d/max_out_d);
sgtitle(sprintf('Test Disturbo: Attenuazione a regime \\approx %.1f dB', att_d_db));

saveas(gcf, fullfile(output_dir, 'reiezione_disturbo.jpg'));


% 2. Rumore di misura n(t) -- MODIFICATO: Durata aumentata e passo ottimizzato
figure(7); clf; set(gcf, 'Color', 'w');

% Passo 1e-6 (sufficiente per 10^5 rad/s) e durata 0.1s (2x tempo assestamento)
tn = 0:1e-6:0.1; 

n_val = zeros(size(tn));
for k = 1:4
    n_val = n_val + 2.5 * sin(1e5 * k * tn);
end

% Risposta al rumore (trasferimento -F)
sim_n = lsim(-F, n_val, tn);

subplot(2,1,1);
plot(tn, n_val, 'r--', 'LineWidth', 0.5); % Linea più sottile per chiarezza
ylabel('Ampiezza'); title('Ingresso: Rumore n(t)');
grid on; xlim([0 max(tn)]);

subplot(2,1,2);
plot(tn, sim_n, 'b', 'LineWidth', 1);
ylabel('Ampiezza'); xlabel('Tempo (s)');
title('Uscita: Residuo y_n(t) (Attenuato)');
grid on; xlim([0 max(tn)]);

% Calcolo attenuazione a regime (ultimi 50% dei campioni)
idx_regime_n = round(length(tn)/2):length(tn);
max_in_n = max(abs(n_val(idx_regime_n)));
max_out_n = max(abs(sim_n(idx_regime_n)));
attenuazione_reale = 20*log10(max_in_n/max_out_n);
sgtitle(sprintf('Test Rumore: Attenuazione a regime \\approx %.1f dB', attenuazione_reale));

saveas(gcf, fullfile(output_dir, 'reiezione_rumore.jpg'));

fprintf('Tutti i grafici aggiornati e corretti (Sforzo Controllo + Rumore esteso).\n');