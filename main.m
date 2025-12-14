clear; close all; clc;

%% 1. Chargement du signal de parole
load('fcno03fz.mat');
x = fcno03fz;
fs = 8000;                     % fréquence d'échantillonnage
x = x / max(abs(x));           % normalisation

%% 2. Ajout d'un bruit blanc gaussien
RSB_dB = 10;                   % rapport signal/bruit en dB
Px = mean(x.^2);
sigma_b = sqrt(Px / (10^(RSB_dB/10)));
b = sigma_b * randn(size(x));
y = x + b;                     % signal bruité

%% 3. Paramètres STFT
Nfft = 1024;
L = 1024;                      % longueur de trame
R = L/2;                       % recouvrement 50%
w = hann(L, 'periodic');       % fenêtre de Hann

%% 4. Estimation du spectre de bruit
nb_trames_bruit = 20;
[B,~,~] = spectrogram(b(1:nb_trames_bruit*R), w, R, Nfft, fs);
Pbb = mean(abs(B).^2, 2);      % densité spectrale moyenne du bruit

%% 5. Calcul STFT (transformée de Fourier à court terme) du signal bruité
[Y, f, t] = spectrogram(y, w, R, Nfft, fs);
Py = abs(Y).^2;

%% 6. Soustraction spectrale
Ps_rehausse = max(Py - Pbb, 0);            % half-wave rectification
S_rehausse_module = sqrt(Ps_rehausse);
S_rehausse = S_rehausse_module .* exp(1j * angle(Y));

%% 7. Sélection d'une trame à observer
num_trame = 60;  % trame choisie pour analyse
start_idx = (num_trame-1)*R + 1;
stop_idx = start_idx + L - 1;

x_frame = x(start_idx:stop_idx) .* w;
y_frame = y(start_idx:stop_idx) .* w;

% trame rehaussée par IFFT du spectre correspondant
S_frame = S_rehausse(:, num_trame);
x_reh_frame = real(ifft(S_frame, Nfft));
x_reh_frame = x_reh_frame(1:L) .* w;

%% 8. AFFICHAGES — Analyse trame par trame
figure('Name','Analyse trame par trame','Color','w');
subplot(2,1,1);
hold on;
plot(x_frame, 'b', 'LineWidth', 1.2);
plot(y_frame, 'r', 'LineWidth', 1.2);
plot(x_reh_frame, 'g', 'LineWidth', 1.2);
title(sprintf('Trame %d : signaux temporels', num_trame));
xlabel('Échantillons'); ylabel('Amplitude');
legend('Original','Bruitée','Rehaussée');
grid on; axis tight;

subplot(2,1,2);
Nfft_local = 1024;
f_axis = linspace(0, fs/2, Nfft_local/2);
Xf = abs(fft(x_frame, Nfft_local));
Yf = abs(fft(y_frame, Nfft_local));
Xrehf = abs(fft(x_reh_frame, Nfft_local));
plot(f_axis, 20*log10(Xf(1:Nfft_local/2)/max(Xf)), 'b', 'LineWidth', 1.2); hold on;
plot(f_axis, 20*log10(Yf(1:Nfft_local/2)/max(Yf)), 'r', 'LineWidth', 1.2);
plot(f_axis, 20*log10(Xrehf(1:Nfft_local/2)/max(Xrehf)), 'g', 'LineWidth', 1.2);
title('Spectres d''amplitude des trois trames');
xlabel('Fréquence (Hz)'); ylabel('Amplitude (dB)');
legend('Original','Bruitée','Rehaussée');
grid on; axis tight;

%% 9. Reconstruction globale par addition-recouvrement
L = length(w);
y_rehausse = zeros(length(y) + L, 1);
poids = zeros(length(y) + L, 1);

for m = 1:length(t)
    frame = real(ifft(S_rehausse(:, m), Nfft));
    frame = frame(1:L) .* w;
    start_idx = (m-1)*R + 1;
    stop_idx = start_idx + L - 1;
    if stop_idx > length(y_rehausse)
        break
    end
    y_rehausse(start_idx:stop_idx) = y_rehausse(start_idx:stop_idx) + frame;
    poids(start_idx:stop_idx) = poids(start_idx:stop_idx) + w.^2;
end
poids(poids < 1e-6) = 1;
y_rehausse = y_rehausse ./ poids;
y_rehausse = y_rehausse(1:length(x));
y_rehausse = y_rehausse / max(abs(y_rehausse));

%% 10. Visualisation comparative : signal original vs rehaussé
figure('Name','Comparaison des signaux original et rehaussé','Color','w');

t = (0:length(x)-1) / fs;

% --- (a) Signaux temporels ---
subplot(2,2,1);
plot(t,x, 'b'); 
title('Signal de parole original');
xlabel('Temps (s)'); ylabel('Amplitude');
grid on; axis tight;

subplot(2,2,2);
plot(t,y_rehausse, 'r');
title('Signal de parole rehaussé');
xlabel('Temps (s)'); ylabel('Amplitude');
grid on; axis tight;

% --- (b) Spectrogrammes ---
subplot(2,2,3);
spectrogram(x, w, R, Nfft, fs, 'yaxis');
title('Spectrogramme du signal original');
xlabel('Temps (s)'); ylabel('Fréquence (kHz)');

subplot(2,2,4);
spectrogram(y_rehausse, w, R, Nfft, fs, 'yaxis');
title('Spectrogramme du signal rehaussé');
xlabel('Temps (s)'); ylabel('Fréquence (kHz)');



%% 11. Évaluation du gain en RSB
RSB_avant = 10*log10(sum(x.^2) / sum((y - x).^2));
RSB_apres = 10*log10(sum(x.^2) / sum((y_rehausse - x).^2));
gain_RSB = RSB_apres - RSB_avant;

fprintf('RSB avant : %.2f dB | RSB après : %.2f dB | Gain : %.2f dB\n', ...
    RSB_avant, RSB_apres, gain_RSB);
