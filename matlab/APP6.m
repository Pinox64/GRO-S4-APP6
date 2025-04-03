clear all;
close all;
[y,Fs] = audioread("note_guitare_LAd.wav");
info = audioinfo("note_guitare_LAd.wav");
duree_son = info.Duration;
N = 160000;

w_hamm = hamming(N);
y_f = y .* w_hamm;
Y = fftshift(fft(y_f));

Fmag = abs(Y);
Fphase = angle(Y);

n = linspace(-Fs/2, Fs/2, N);

%% Trouver harmoniques
f0 = 466;
harmoniques = f0*(1:32);

index_harmo = zeros(1, length(harmoniques));
for k = 1:length(harmoniques)
    % Chercher l'index le plus proche
    [~, idx] = min(abs(n - harmoniques(k)));
    
    % Chercher localement le maximum dans une petite fenêtre
    range = max(1, idx-1000):min(N, idx+1000); % Éviter les dépassements
    [~, local_max] = max(Fmag(range));
    
    % Mettre à jour l'index avec la position du vrai pic
    index_harmo(k) = range(local_max);
end

figure(1);
plot(n, 20*log10(Fmag));
hold on;
scatter(n(index_harmo), 20*log10(Fmag(index_harmo)), 'ro', 'filled'); 
xlabel('Fréquence Hz');
ylabel('Magnitude (dB)');
title('Détection des harmoniques');
grid on;
legend('FFT', 'Harmoniques détectées');
hold off;
figure(2);
plot(n ,Fphase);
xlabel('Fréquence Hz');
ylabel('Phase');

%% Passe-bas RIF

%Ordre du filtre
Fc = pi/1000;
N = 1000;
m = N*Fc/Fs;
K = 2*m+1;

% Génération de la réponse impulsionnelle
k = -N/2:N/2-1; % Indices centrés autour de 0
h = zeros(size(k)); % Initialisation du filtre

% Calcul des coefficients
for i = 1:length(k)
    if k(i) == 0
        h(i) = K / N;
    else
        h(i) = (1/N) * (sin(pi * k(i) * K / N) / sin(pi * k(i) / N));
    end
end

h = hamming(N)'.*h;
y_abs = abs(y);

y_filtered = conv(y_abs, h, 'same');

figure(3);
plot(y_abs, 'b'); hold on;
plot(y_filtered, 'r'); hold off;
xlabel('Échantillon');
ylabel('Magnitude');
legend('Avant filtrage', 'Après filtrage');
title('Signal et enveloppe temporelle');
grid on;

%%Recréer LA#

t = (0:length(y_filtered)-1) / Fs;

sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthLA = sum_sinuses' .* y_filtered;
synthLA = synthLA / max(abs(synthLA));

%figure(4);
%plot(synthLA);
%hold on;
%plot(y);
%hold off;

%Synth LA#
%sound(synthLA, Fs);

%%Recréer SOL
t = (0:length(y_filtered)-1) / Fs;

SOL_harmoniques = harmoniques*0.891;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*SOL_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthSOL = sum_sinuses' .* y_filtered;
synthSOL = synthSOL / max(abs(synthSOL));

%T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',4);
%start(T);
%wait(T);
%sound(synthSOL, Fs);

%%Recréer MI
t = (0:length(y_filtered)-1) / Fs;

MI_harmoniques = harmoniques*0.749;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*MI_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthMI = sum_sinuses' .* y_filtered;
synthMI = synthMI / max(abs(synthMI));

%T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',4);
%start(T);
%wait(T);
%sound(synthMI, Fs);

%%Recréer FA
t = (0:length(y_filtered)-1) / Fs;

FA_harmoniques = harmoniques*0.794;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*FA_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthFA = sum_sinuses' .* y_filtered;
synthFA = synthFA / max(abs(synthFA));

%T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',4);
%start(T);
%wait(T);
%sound(synthFA, Fs);

%%Recréer RE
t = (0:length(y_filtered)-1) / Fs;

RE_harmoniques = harmoniques*0.667;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*RE_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthRE = sum_sinuses' .* y_filtered;
synthRE = synthRE / max(abs(synthRE));

%T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',4);
%start(T);
%wait(T);
%sound(synthRE, Fs);

%% Jouer mélodie complète

T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',0.25);
T2 = timer('TimerFcn',@(~,~)disp(''),'StartDelay',2);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthMI, Fs);
start(T2);
wait(T2);
stop(T2);

sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthRE, Fs);

