function mini_interface_soustraction_spectrale

    fig = uifigure('Name', 'Rehaussement de parole', ...
                   'Position', [100 100 600 400]);

    % Bouton pour charger le signal
    uibutton(fig, 'Text', 'Charger un signal', ...
        'Position', [50 340 150 30], ...
        'ButtonPushedFcn', @(btn,event) charger_signal());

    % Slider RSB
    uilabel(fig, 'Text', 'RSB (dB)', 'Position', [250 340 60 30]);
    sld = uislider(fig, 'Position', [320 355 200 3], ...
                   'Limits', [0 30], 'Value', 10);

    % Bouton traitement
    uibutton(fig, 'Text', 'Appliquer traitement', ...
        'Position', [50 280 150 30], ...
        'ButtonPushedFcn', @(btn,event) traiter_signal(sld));

    % Axes spectrogramme
    ax = uiaxes(fig, 'Position', [50 50 500 200]);
    title(ax, 'Spectrogramme du signal');

    % ============================== FONCTIONS ==============================
    function charger_signal()
        [f,p] = uigetfile({'*.mat;*.wav'});
        if f ~= 0
            assignin('base','signal_path',fullfile(p,f));
            uialert(fig,'Signal chargé avec succès','OK');
        end
    end

    function traiter_signal(slider)
        try
            if ~evalin('base','exist(''signal_path'',''var'')')
                uialert(fig,'Veuillez charger un signal.','Erreur');
                return;
            end

            RSB = slider.Value;
            signal_path = evalin('base','signal_path');
            [~,~,ext] = fileparts(signal_path);

            % Chargement
            if strcmp(ext,'.mat')
                s = load(signal_path);
                x = struct2array(s);
                fs = 8000;
            else
                [x,fs] = audioread(signal_path);
            end

            x = x(:);
            x = x / max(abs(x));

            % Bruitage
            Px = mean(x.^2);
            sigma_b = sqrt(Px / (10^(RSB/10)));
            y = x + sigma_b*randn(size(x));

            % ====== PARAMÈTRES STFT ======
            Nfft = 1024;
            w = hann(Nfft,'periodic');
            R = Nfft/2;

            % ========= ESTIMATION BRUIT =========
            nb_trames_bruit = 20;
            Lb = nb_trames_bruit * R;

            % FORÇAGE : minimum = 2 trames complètes
            if Lb < 2*Nfft
                Lb = 2*Nfft;
            end
            if Lb > length(y)
                Lb = length(y);
            end

            [B,~,~] = spectrogram(y(1:Lb), w, R, Nfft, fs);

            % GARANTIR QUE B A 1024 LIGNES
            if size(B,1) ~= Nfft
                error('STFT du bruit incorrecte (%d lignes).', size(B,1));
            end

            Pbb = mean(abs(B).^2,2);

            % ========= STFT SIGNAL BRUITÉ =========
            [Y,~,~] = spectrogram(y, w, R, Nfft, fs);
            Py = abs(Y).^2;

            % ========= SOUSTRACTION =========
            alpha = 0.98;
            floor_gain = 0.02;
            Ps = max(Py - alpha*Pbb, floor_gain*Pbb);
            S = sqrt(Ps).*exp(1j*angle(Y));

            % ========= RECONSTRUCTION =========
            y_out = istft(S, fs, ...
                          'Window', w, ...
                          'OverlapLength', R, ...
                          'FFTLength', Nfft);

            y_out = y_out / max(abs(y_out));

            % ========= AFFICHAGE =========
            spectrogram(y_out,w,R,Nfft,fs,'yaxis','Parent',ax);
            title(ax, sprintf('Signal rehaussé (RSB = %.1f dB)',RSB));

            % ========= LECTURE =========
            soundsc(y_out, fs);

        catch ME
            uialert(fig, ME.message, 'Erreur STFT');
        end
    end
end
