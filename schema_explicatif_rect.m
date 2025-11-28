%% Overlap–Add avec fenêtre RECTANGULAIRE (2 trames)
clear; close all; clc;

L = 40;
R = L/2;
N = L + R;
n = 0:N-1;

x = sin(2*pi*n/40);
w = rectwin(L)';

%% Trames
x1 = x(1:L).*w;
x2 = x(R+1:R+L).*w;

%% Reconstruction
x_recon = zeros(1,N);
weight  = zeros(1,N);

x_recon(1:L) = x1;
x_recon(R+1:R+L) = x_recon(R+1:R+L) + x2;

weight(1:L) = w;
weight(R+1:R+L) = weight(R+1:R+L) + w;

x_recon_norm = x_recon ./ weight;

%% FIGURE
figure('Color','w','Position',[80 80 900 750]);

subplot(3,1,1); hold on;
plot(n,x,'k','LineWidth',1.3);

area(0:L-1, x1,'FaceColor',[0.3 0.6 1],'FaceAlpha',0.4,'EdgeColor','none');
area(R:R+L-1, x2,'FaceColor',[1 0.6 0.3],'FaceAlpha',0.4,'EdgeColor','none');

w_full = [0, w, 0];

% --- Fenêtres affichées en transparence ---
patch([-1:L], w_full, [1 0 0], ...
    'FaceAlpha', 0.25, 'EdgeColor', 'r', 'LineWidth', 1.4);

patch([R-1:R+L], w_full, [0 0.3 1], ...
    'FaceAlpha', 0.25, 'EdgeColor', [0 0.3 1], 'LineWidth', 1.4);

title('(a) Fenêtrage RECTANGULAIRE');
xlabel('n'); ylabel('Amplitude');
grid on; axis tight;


subplot(3,1,2); hold on;
yyaxis left
plot(n, x_recon,'m','LineWidth',1.3);
ylabel('Somme des trames');

yyaxis right
plot(n, weight,'g--','LineWidth',1.3);
ylabel('\Sigma w[n]');

title('(b) Somme des trames (Rectwin)');
xlabel('n');
grid on; axis tight;


subplot(3,1,3); hold on;
plot(n,x,'k','LineWidth',1.3);
plot(n,x_recon_norm,'r--','LineWidth',1.3);
plot(n,x-x_recon_norm,'b:','LineWidth',1.2);

title('(c) Reconstruction (Rectwin)');
xlabel('n'); ylabel('Amplitude');
legend('x[n]','x_{recon}','Erreur');
grid on; axis tight;

sgtitle('Procédure Overlap–Add — Fenêtre RECTANGULAIRE');
