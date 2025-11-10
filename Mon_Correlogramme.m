function [P_corr, f] = Mon_Correlogramme(x, Nfft, Fe)
    % Corrélogramme : estimation de la DSP à partir de l'autocorrélation
    [R, ~] = xcorr(x, 'biased');     % autocorrélation estimée
    P_corr = abs(fftshift(fft(R, Nfft))) ;  % FFT de Rxx
    f = linspace(-Fe/2, Fe/2, Nfft);
end
