function [P_smooth, f] = Mon_Daniell(x, Nfft, Fe, L)
    % Périodogramme de Daniell : lissage fréquentiel du périodogramme
    %   L    : taille du lissage (nombre de points de moyenne)
    
    [P, f] = Mon_Periodogramme(x, Nfft, Fe);
    
    % Lissage fréquentiel 
    win = ones(1, L) / L;
    P_smooth = conv(P, win, 'same');
end
