clear all
close all
load("data_basson.mat")
%%Recréer un LA#
t = (0:length(y)-1) / fe;

sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthLA = sum_sinuses' .* y_enveloppe;
synthLA = synthLA / max(abs(synthLA));

%sound(synthLA, fe);
audiowrite("basson_generated.wav",synthLA, fe);

% Affichage de la comparaison des deux signaux (original vs généré)
figure(29);
clf
fftSynth = abs(fftshift(fft(synthLA,N)));
subplot(3,1,3);
plot(n,20*log(fftSynth))
title("Fréquences du signal basson synthétisé");
xlabel("fréquence (Hz)");
ylabel("Amplitude (dB)");
subplot(3,1,2);
fftFiltr = abs(fftshift(fft(y, N)));
plot(n,20*log(fftFiltr));
title("Fréquences du signal basson filtré");
xlabel("fréquence (Hz)");
ylabel("Amplitude (dB)");
subplot(3,1,1)
[x,fe] = audioread("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs.wav");
fftOriginal = abs(fftshift(fft(x,N)));
plot(n,20*log(fftOriginal));
title("Fréquences du signal basson original");
xlabel("fréquence (Hz)");
ylabel("Amplitude (dB)");