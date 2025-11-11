%% --- PRELIMINAIRE 4 : Addition - Recouvrement (Overlap-Add) ---
clear; close all; clc;

%% 1. Création d'un signal synthétique
fs = 1600;                    % fréquence d'échantillonnage
t = (0:1/fs:1-1/fs)';         % vecteur temps
x = sin(2*pi*50*t) + 0.5*sin(2*pi*120*t); %signal synthétique multi-fréquence


%% 2. Paramètres de tramage
L = 15;           % longueur de trame
R = L/2;            % pas de trame (50% de recouvrement)
w = hann(L, 'periodic');  % fenêtre de Hann

%% 3. Découpage du signal en trames fenêtrées
N = length(x);
M = floor((N - L)/R) + 1;   % nombre de trames
frames = zeros(L, M);

for m = 1:M
    start_idx = (m-1)*R + 1;
    frames(:,m) = x(start_idx:start_idx+L-1) .* w;
end

%% 4. Reconstruction par addition-recouvrement
x_recon = zeros(N,1);
weight  = zeros(N,1);

for m = 1:M
    start_idx = (m-1)*R + 1;
    x_recon(start_idx:start_idx+L-1) = x_recon(start_idx:start_idx+L-1) + frames(:,m);
    weight(start_idx:start_idx+L-1)  = weight(start_idx:start_idx+L-1) + w;
end

% éviter division par zéro
weight(weight==0) = 1;
x_recon = x_recon ./ weight;

%% 5. Vérification de la reconstruction
erreur = x - x_recon;
EQM = mean(erreur.^2);

fprintf('Erreur quadratique moyenne : %.2e\n', EQM);
fprintf('Amplitude maximale de l''erreur : %.2e\n', max(abs(erreur)));


%% 6. Affichages
figure;
subplot(3,1,1);
plot(x); title('Signal original');

subplot(3,1,2);
plot (weight, 'LineWidth', 1.2);
title('Somme des fenêtres (poids de recouvrement)');
xlabel('Échantillons n'); ylabel('Somme des pondérations');
grid on;

subplot(3,1,3);
plot(x_recon, '--r'); title('Signal reconstruit');

figure;
plot(erreur, 'k'); hold on;
yline(0, 'r--');
xlim([L, N-L]); % zoom sur la zone centrale
title('Erreur de reconstruction (zone centrale)');
grid on;

