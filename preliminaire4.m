%% --- PRELIMINAIRE 4 : Addition - Recouvrement (Overlap-Add) ---
clear; close all; clc;

%% 1. Création d'un signal synthétique
fs = 1600;                          % fréquence d'échantillonnage
t  = (0:1/fs:1-1/fs)';              % vecteur temps
x  = sin(2*pi*50*t) + 0.5*sin(2*pi*120*t);   % signal multi-fréquence
%x = x .* (1 + 0.4*randn(size(x)));  % rendre le signal non stationnaire

N  = length(x);

%% ================== CAS 1 : FENETRE DE HANN ==================
L_hann = 32;                        % longueur de trame (pair)
R_hann = L_hann/8;                  % pas de trame (50 %)
w_hann = hann(L_hann);              % fenêtre de Hann (NON périodique)

% --- Découpage ---
M_hann    = floor((N - L_hann)/R_hann) + 1;
frames_h  = zeros(L_hann, M_hann);

for m = 1:M_hann
    start_idx = (m-1)*R_hann + 1;
    frames_h(:,m) = x(start_idx:start_idx+L_hann-1) .* w_hann;
end

% --- Reconstruction OLA ---
x_recon_h    = zeros(N,1);
weight_hann  = zeros(N,1);

for m = 1:M_hann
    start_idx = (m-1)*R_hann + 1;
    x_recon_h(start_idx:start_idx+L_hann-1) = x_recon_h(start_idx:start_idx+L_hann-1) + frames_h(:,m);
    weight_hann(start_idx:start_idx+L_hann-1)  = weight_hann(start_idx:start_idx+L_hann-1) + w_hann;
end

weight_hann(weight_hann==0) = 1;
x_recon_h = x_recon_h ./ weight_hann;

% --- Erreur ---
erreur_hann = x - x_recon_h;
EQM_hann    = mean(erreur_hann.^2);

fprintf('--- Fenêtre HANN ---\n');
fprintf('EQM : %.2e\n', EQM_hann);
fprintf('Amplitude max erreur : %.2e\n\n', max(abs(erreur_hann)));

%% ================= CAS 2 : FENETRE RECTANGULAIRE =================
L_rect = 32;                        % même longueur pour comparer
R_rect = L_rect/8
w_rect = ones(L_rect,1);            % fenêtre RECTANGULAIRE

% --- Découpage ---
M_rect   = floor((N - L_rect)/R_rect) + 1;
frames_r = zeros(L_rect, M_rect);

for m = 1:M_rect
    start_idx = (m-1)*R_rect + 1;
    frames_r(:,m) = x(start_idx:start_idx+L_rect-1) .* w_rect;
end

% --- Reconstruction OLA ---
x_recon_r   = zeros(N,1);
weight_rect = zeros(N,1);

for m = 1:M_rect
    start_idx = (m-1)*R_rect + 1;
    x_recon_r(start_idx:start_idx+L_rect-1) = x_recon_r(start_idx:start_idx+L_rect-1) + frames_r(:,m);
    weight_rect(start_idx:start_idx+L_rect-1)  = weight_rect(start_idx:start_idx+L_rect-1) + w_rect;
end

x_recon_h_no_norm = x_recon_h .* 0 + x_recon_h; % copie
x_recon_r_no_norm = x_recon_r .* 0 + x_recon_r; % copie
erreur_hann_no_norm = x - x_recon_h_no_norm;
erreur_rect_no_norm = x - x_recon_r_no_norm;


weight_rect(weight_rect==0) = 1;
x_recon_r = x_recon_r ./ weight_rect;



% --- Erreur ---
erreur_rect = x - x_recon_r;
EQM_rect    = mean(erreur_rect.^2);

fprintf('--- Fenêtre RECTANGULAIRE ---\n');
fprintf('EQM : %.2e\n', EQM_rect);
fprintf('Amplitude max erreur : %.2e\n', max(abs(erreur_rect)));

% --- Hann : signal / poids / recon ---
figure;
subplot(3,1,1);
plot(x); title('Signal original');
grid on;

subplot(3,1,2);
plot(weight_hann,'LineWidth',1.2);
title('Somme des fenêtres (Hann)');
xlabel('Échantillons'); ylabel('Poids');
grid on;

subplot(3,1,3);
plot(x_recon_h,'--r');
title('Signal reconstruit (fenêtre de Hann)');
grid on;

% --- Rect : signal / poids / recon ---
figure;
subplot(3,1,1);
plot(x); title('Signal original');
grid on;

subplot(3,1,2);
plot(weight_rect,'LineWidth',1.2);
title('Somme des fenêtres (rectangulaire)');
xlabel('Échantillons'); ylabel('Poids');
grid on;

subplot(3,1,3);
plot(x_recon_r,'--r');
title('Signal reconstruit (fenêtre rectangulaire)');
grid on;

% --- Erreur Hann seule ---
figure;
plot(erreur_hann,'k'); hold on;
yline(0,'r--');
xlim([L_hann, N-L_hann]);
title('Erreur de reconstruction (Hann) - zone centrale');
grid on;

% --- Erreur Rect seule ---
figure;
plot(erreur_rect,'k'); hold on;
yline(0,'r--');
xlim([L_rect, N-L_rect]);
title('Erreur de reconstruction (Rectangulaire) - zone centrale');
grid on;




%% === COMPARAISON DES SPECTRES ===
Nfft = 4096;

Xh = abs(fft(x_recon_h, Nfft));
Xr = abs(fft(x_recon_r, Nfft));
Xo = abs(fft(x, Nfft));

f = linspace(0, fs/2, Nfft/2);

figure;
plot(f,20*log10(Xo(1:Nfft/2)),'k','LineWidth',1.3); hold on;
plot(f,20*log10(Xh(1:Nfft/2)),'b--','LineWidth',1.2);
plot(f,20*log10(Xr(1:Nfft/2)),'r--','LineWidth',1.2);
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
legend('Original','Reconstruit Hann','Reconstruit Rect','Location','best');
title('Comparaison spectrale des signaux reconstruits');
grid on;

figure;
subplot(2,1,1);
plot(erreur_hann_no_norm);
title('Erreur SANS normalisation — Hann');
grid on;

subplot(2,1,2);
plot(erreur_rect_no_norm);
title('Erreur SANS normalisation — Rectwin');
grid on;
