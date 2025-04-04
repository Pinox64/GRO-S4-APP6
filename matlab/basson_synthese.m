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

sound(synthLA, fe);
audiowrite("basson_generated.wav",synthLA, fe);

% Affichage de la comparaison des deux signaux (original vs généré)
figure(29);
clf
fftSynth = abs(fftshift(fft(synthLA,N)));
subplot(2,1,1);
plot(n,20*log(fftSynth))
title("Généré");
subplot(2,1,2);
fftOriginal = abs(fftshift(fft(y)));
plot(n,20*log(fftOriginal));
title("original");