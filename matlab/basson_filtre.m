clear all
close all
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
w0 = 2*pi*fc_nf/fe; %Fréqence centrale (en rad/échantillons) de la bande à élminier (
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

Nfft = length(x)*20;
n_ = -N/2 : N/2-1; 
h_ = hamming(N)'.*h(n_);   % Reponse impulsionnelle du filtre coupe bande avec effet de gibs supprimé (fenétrage)
hb_ = h_bas2(n_); 
h_fft = fftshift(fft(h_,Nfft));
h_mag = 20*log(abs(h_fft));
h_phase = angle(h_fft);

hb_fft = fftshift(fft(hb_,Nfft));
hb_mag = 20*log(abs(hb_fft));
hb_phase = angle(hb_fft);

f_ = linspace(-fe/2, fe/2, Nfft);   % axe des fréquences en Hz
%Vérification du filtre en l'appliquant a une sinus de 1000 Hz
t = 0:4*2*pi/w0;
sin1000 = sin(w0*t);
sin1000_filtered = conv(sin1000, h_, "same");
sin1000_2filtered = conv(sin1000_filtered,h_,"same");
sin1000_3filtered = conv(sin1000_2filtered,h_,"same");
sin1000_4filtered = conv(sin1000_3filtered,h_,"same");
sin1000_5filtered = conv(sin1000_4filtered,h_,"same");
figure('Name','Réponse du filtre coupe bande a une sinus de 1000 Hz')
plot(t, sin1000);
hold on
plot(t,sin1000_filtered);
plot(t,sin1000_2filtered);
plot(t,sin1000_3filtered);
plot(t,sin1000_4filtered);
plot(t,sin1000_5filtered);
legend("Entrée", "Signal filtré 1 fois", "Signal filtré 2 fois", "Signal filtré 3 fois", "Signal filtré 4 fois", "Signal filtré 5 fois")
title("Réponse du filtre coupe bande à un signal d'entrée sinusoidal de 1000 Hx")
ylabel("amplitude");
xlabel("échantillons")

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
ylabel("Gain(db)");
xlabel("Fréquence (Hz)");
title("Réponse fréquentielle du filtre RIF coupe bande");
xlim([-2000 2000]);  % Zoom de 900 à 1100 Hz
subplot(2,1,2);
plot(f_,h_phase);
ylabel("Phase (rad)");
xlabel("Fréquence (Hz)");
xlim([-2000 2000]);  % Zoom de 900 à 1100 Hz
%title("Phase");
%legend("Coupe bande", "Passe bas");

% Affichage des poles et des zero du coupe band
figure ("Name", "Pôles et zéros du filtre RIF coupe bande");
zplane(h_,1);
title("Pôles et zéros du filtre RIF coupe bande");

figure ("Name","Coupe bande");
clf
plot(n_,h_);
title("Réponse impulsionnelle temporelle du filtre coupe bande")
xlabel("Échantillons")
ylabel("Magnitude")
%h_trunc = h_(end/2, end)
h_mogus = hamming(N)'.*h_;
filteredsinSound = conv(conv(x,h_mogus, "same"),h_mogus,"same"); % On applique le filtre 2 fois ;)
audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_NoSin.wav",filteredsinSound,fe)

% Graph entrée vs sortie 1000 Hz
fftIn = abs(fftshift(fft(x,Nfft)));
fftOut = abs(fftshift(fft(filteredsinSound,Nfft)));
figure("Name","Entrée vs Sortie filtre 1000 Hz")
clf
subplot(2,1,1);
plot(f_, 20*log(fftIn));
title("Spectres du signal d'entrée");
xlim([-2000 2000])
xlabel("Fréquence (Hz)");
ylabel("Magnitude (dB)");
subplot(2,1,2);
plot(f_,20*log(fftOut));
title("Spectres du signal filtré par le coupe bande");
xlim([-2000 2000])
xlabel("Fréquence (Hz)");
ylabel("Magnitude (dB)");

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
 fullyFiltered = filtfilt(b,a,filteredsinSound);
 audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_hf_couper.wav",fullyFiltered,fe);
 %figure (1);
 %t = 1:length(x);
 %plot(t,x);
 %hold on;
 %plot(t,fullyFiltered);
 %legend("in","out");

 % Affichage des poles et des zero du filtre valant du beurre
 figure ("Name", "Pôles et zéros du filtre RII (butterworth)");
 zplane(b,a);
 title("Pôles et zéros du filtre RII (butterworth)");
 figure("Name","Filtre butter")
 freqz(b,a)


% Graph entrée vs sortie butter
hammingSurFiltered = fullyFiltered.*hamming(length(fullyFiltered));
fftOutButter = abs(fftshift(fft(fullyFiltered,Nfft)));
figure("Name","Entrée vs Sortie filtre Butterworth")
clf
subplot(2,1,2);
plot(f_, 20*log(fftOutButter));
title("Spectres du signal filtré par le butterworth");
%xlim([-2000 2000])
xlabel("Fréquence (Hz)");
ylabel("Magnitude (dB)");
subplot(2,1,1);
plot(f_,20*log(fftOut));
title("Spectres du signal d'entrée");
%xlim([-2000 2000])
xlabel("Fréquence (Hz)");
ylabel("Magnitude (dB)");

%% Rééchantillonnement du signal Basson
y_downsampled = downsample(fullyFiltered,2);
fe2 = fe/2
N = length(y_downsampled);
n = linspace(-fe2/2, fe2/2, N);

%sound(y_downsampled, fe2)
audiowrite("note_basson_filtered_and_downsampled.wav", y_downsampled, fe2)