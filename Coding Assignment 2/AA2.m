close all; clear

N = 220;
cutoff = 0.325;

h_blackman = fir1(N-1, cutoff, blackman(N));
[H_blackman, w] = freqz(h_blackman, 1, 8000);

figure;
plot(w/pi, abs(H_blackman), 'r', 'DisplayName', 'Blackman');
hold on;
plot(w/pi, double(abs(w) <= 1.02), 'k--', 'DisplayName', 'Ideal LPF');
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude');
legend;
title('FIR Filter vs Ideal Filter');
grid on;
