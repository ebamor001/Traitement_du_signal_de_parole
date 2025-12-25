clear; close all; clc;

%% Chargement du signal
load fcno03fz.mat
x = fcno03fz(:);
fs = 8000;
x = x / max(abs(x));   

%% Paramètres STFT
L = 1024;
R = L/2;
Nfft = 1024;
w = hann(L, "periodic");

%% Paramètres d'évaluation
RSB_tests = [5 10 15];     % RSB à tester
Nreal = 20;                % nombre de réalisations
results = zeros(length(RSB_tests), 3); % [RSB_avant, RSB_apres, Gain]

for k = 1:length(RSB_tests)
    RSBdB = RSB_tests(k);
    gains = zeros(Nreal,1);

    for r = 1:Nreal
        
        %% --- Génération du bruit ---
        Px = mean(x.^2);
        sigma = sqrt(Px / 10^(RSBdB/10));
        b = sigma * randn(size(x));
        y = x + b;

        %% --- Estimation du bruit ---
        nb_trames_bruit = 20;
        [B,~,~] = spectrogram(b(1:nb_trames_bruit*R), w, R, Nfft, fs);
        Pbb = mean(abs(B).^2, 2);

        %% --- STFT du signal bruité ---
        [Y,~,~] = spectrogram(y, w, R, Nfft, fs);
        Py = abs(Y).^2;

        %% --- Soustraction spectrale ---
        Ps_est = max(Py - Pbb, 0);
        S_module = sqrt(Ps_est);
        S_hat = S_module .* exp(1j*angle(Y));

        %% --- Reconstruction par OLA ---
        y_hat = zeros(length(y) + L, 1);
        poids = zeros(length(y) + L, 1);

        for m = 1:size(S_hat,2)
            frame = real(ifft(S_hat(:,m), Nfft));
            frame = frame(1:L).*w;
            idx = (m-1)*R+1 : (m-1)*R+L;
            if idx(end) > length(y_hat), break; end
            y_hat(idx) = y_hat(idx) + frame;
            poids(idx) = poids(idx) + w.^2;
        end
        
        poids(poids<1e-6) = 1;
        y_hat = y_hat ./ poids;
        y_hat = y_hat(1:length(x));

        %% --- Calcul RSB ---
        RSB_avant  = 10*log10( sum(x.^2) / sum((y      - x).^2) );
        RSB_apres  = 10*log10( sum(x.^2) / sum((y_hat  - x).^2) );
        gains(r)   = RSB_apres - RSB_avant;
    end

    %% --- Résultats moyens ---
    results(k,:) = [RSB_tests(k), mean(RSB_apres), mean(gains)];
end

disp("=== Tableau RSB (moyenne sur "+Nreal+" réalisations) ===")
disp("RSB_avant | RSB_apres | Gain_RSB")
disp(results)
