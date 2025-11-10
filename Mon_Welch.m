function [P_welch, f] = Mon_Welch(x, Nfft, Fe,overlap)
    hop = round(Nfft * (1 - overlap));  % pas entre segments
    
    N = length(x);
    nb_segments = floor((N - Nfft) / hop) + 1;
    P_welch = zeros(1, Nfft);
    
    for i = 1:nb_segments
        % extraire le segment
        idx = (i-1)*hop + 1 : (i-1)*hop + Nfft;
        seg = x(idx) ;
    
        % calcul du périodogramme
        [P, ~] = Mon_Periodogramme(seg, Nfft, Fe);
    
        % accumulation
        P_welch = P_welch + P;
    end
    
    % moyenne
    P_welch = P_welch / nb_segments;
    
    % centrage et axe fréquentiel
    P_welch = fftshift(P_welch);
    f = linspace(-Fe/2, Fe/2, Nfft);
end
