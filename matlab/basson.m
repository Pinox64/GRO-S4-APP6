clear all

%But: filtrer le bruit et la sinus.
%Pour filtrer le bruit on veut faire un RII un passe bas pour eliminer +15khz
%Pour filtrer la sinus, on veut faire un Coupe bande butterworth pour
%eliminer tout ce qui est entre 980hz à 1020hz
[x,fe] = audioread("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs.wav");
 %% Filtre coupe bande

fc_lp = 40;  % Fréquence de coupure du filtre
N = 1024;   %Ordre du filtre (Donné dans l'énoncé)
m = N*fc_lp/fe
K = m*2+1;
fc_nf = 1000;
w0 = 2*pi*fc_nf/fe; %Fréqence centrale (en Hz) de la bande à élminier (
% A small local function that never computes 0/0 for x=0
function val = h_bas(x, K, N)
    if x == 0
        val = K/N;
    else
        val = (1/N) * sin(pi*K*x/N) / sin(pi*x/N);
    end
end
h_bas2 = @(n) arrayfun(@(x) h_bas(x, K, N), n);
function val = h_bande(x, K, N,w0)
    delta = double(x==0); %Dirac discret
    val = delta - 2*h_bas(x,K,N)*cos(w0*x);
end
% Then define an anonymous "vectorized" function using arrayfun:
h = @(n) arrayfun(@(x) h_bande(x, K, N,w0), n);

Nfft = N;
n_ = -N/2 : N/2-1; 
h_ = h(n_);   % No NaNs, h_ is a vector
hb_ = h_bas2(n_); 
h_fft = fftshift(fft(h_,Nfft));
h_mag = 20*log(abs(h_fft));
h_phase = angle(h_fft);

hb_fft = fftshift(fft(hb_,Nfft));
hb_mag = 20*log(abs(hb_fft));
hb_phase = angle(hb_fft);

f_ = linspace(-fe/2, fe/2, Nfft);   % axe des fréquences en Hz

figure("Name","reponse freq du passe bas")
clf
%freqz(h_,1024);
%hold on;
%freqz(h_b,1024);
subplot(2,1,1);
plot(f_,hb_mag);
title("mag");
subplot(2,1,2);
plot(f_,hb_phase);
title("phase");
xlabel("Frequence");
ylabel("Amplitude (db)");
%legend("Coupe bande", "Passe bas");

figure("Name","reponse freq du coupe bande")
clf
%freqz(h_,1024);
%hold on;
%freqz(h_b,1024);
subplot(2,1,1);
plot(f_,h_mag);
title("mag");
subplot(2,1,2);
plot(f_,h_phase);
title("phase");
legend("Coupe bande", "Passe bas");

figure ("Name","Coupe bande");
clf
plot(n_,h_);
title("Réponse a une impulsion du coupe bande")
%h_trunc = h_(end/2, end)
h_mogus = hamming(Nfft)'.*h_;
h_reconstructed = real(ifft(ifftshift(h_)));
filteredsinSound = conv(x,h_mogus, "same");
filteredsinSound = filteredsinSound(length(h_mogus):length(filteredsinSound)-length(h_mogus));
audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_NoSin.wav",filteredsinSound,fe)

%% Filtre RII (butterworth) pour filtrer ce qui est plus haut que 15KHz

%fc = 12000;
fp = 8000; % Fréquence de passage (on coupe moins de 1db)
fs = 15000; % Fréquence qu'on veut complètement couper
ny = fe/2; % Fréquence de nyquist du signal
Wp = fp/ny; % Fréquence de passage normalisée selon nyquist de fe
Ws = fs/ny; % Fréquence de stop normalisée selon nyquist de fe

Rp = 1; % Ripple minimal dans la fréquence de passage (en db)
Rs = 40; % Atténuation minimale à la fréquence de stop (en db)

[N, Wn] = buttord(Wp,Ws,Rp,Rs)
[b,a] = butter(N,Wn,"low");
 fullyFiltered = filter(b,a,filteredsinSound);
 audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_hf_couper.wav",fullyFiltered,fe);
 %figure (1);
 %t = 1:length(x);
 %plot(t,x);
 %hold on;
 %plot(t,fullyFiltered);
 %legend("in","out");
y = fullyFiltered;
w_hamm = hamming(length(y));
y_f = y .* w_hamm;      %Appliquer un fenétrage hamming sur le signal d'entrée
Y = fftshift(fft(y_f));


%% Rééchantillonnement du signal Basson
y_downsampled = downsample(fullyFiltered,2);
fe2 = fe/2
N = length(y_downsampled);
n = linspace(-fe2/2, fe2/2, N);

sound(y_downsampled, fe2)

%% Trouver les harmoniques du signal Basson et leur magnitude
F = fftshift(fft(y_downsampled,N));
Fmag = abs(F);
Fphase = angle(F);

f0 = 466/2;                                               % Fréquence fondamentale du son d'origine (La#)
harmoniques = f0*(1:32);                                % Fréquence de chaque harmonique théorique de la fondamentale

index_harmo = zeros(1, length(harmoniques));            %Index des maximums de la fft qui correspondent aux harmoniques théoriques
for k = 1:length(harmoniques)
    % Chercher l'index fréquentiel le plus proche de l'harmonique théorique
    [~, idx] = min(abs(n - harmoniques(k)));
    
    % Chercher localement le maximum dans une petite fenêtre
    range = max(1, idx-1000):min(N, idx+1000);          % Fenétrage pour Éviter les dépassements
    [~, local_max] = max(Fmag(range));                  % Trouver l'index fréquentiel dont la magnitude est la plus haute
    
    % Mettre à jour l'index avec la position du vrai pic
    index_harmo(k) = range(local_max);
end

figure(21);
clf
plot(n, 20*log(Fmag));
hold on;
scatter(n(index_harmo), 20*log(Fmag(index_harmo)), 'ro', 'filled'); 
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

%% Passe-bas RIF (génération de l'enveloppe)

%Ordre du filtre
Fc = pi/1000;           % Fréquence de coupure
N_lp = 1000;               % Ordre du filtre
m = N_lp*Fc/fe2;
K = 2*m+1;

% Génération de la réponse impulsionnelle
k = -N_lp/2:N_lp/2-1; % Indices centrés autour de 0
h = zeros(size(k)); % Initialisation du filtre

% Calcul des coefficients de la réponse impulsionnelle
for i = 1:length(k)
    if k(i) == 0
        h(i) = K / N_lp;
    else
        h(i) = (1/N_lp) * (sin(pi * k(i) * K / N_lp) / sin(pi * k(i) / N_lp));
    end
end

h = hamming(N_lp)'.*h;     % On applique un filtre hamming afin d'éviter l'effet de fuite
y_abs = abs(y);         % On met le signal d'entrée en valeur absolue car l'envloppe temporelle est tjrs au dessus de zéro

y_filtered = conv(y_abs, h, 'same');    % On applique le filtre passe bas au signal d'entrée

% Affichage du résultat de l'enveloppe temp
figure(23);
clf
plot(y_abs, 'b'); hold on;
plot(y_filtered, 'r'); hold off;
xlabel('Échantillon');
ylabel('Magnitude');
legend('Avant filtrage', 'Après filtrage');
title('Signal et enveloppe temporelle');
grid on;

%%Recréer un LA#
t = (0:length(y_filtered)-1) / fe2;

sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthLA = sum_sinuses' .* y_filtered;
synthLA = synthLA / max(abs(synthLA));

sound(synthLA, fe2);
audiowrite("basson_generated.wav",synthLA, fe2);

% Affichage de la comparaison des deux signaux (original vs généré)
figure(29);
clf
fftSynth = abs(fftshift(fft(synthLA,N)));
subplot(2,1,1);
plot(n,20*log(fftSynth))
title("Généré");
subplot(2,1,2);
fftOriginal = abs(fftshift(fft(y_downsampled)));
plot(n,20*log(fftOriginal));
title("original");