clear; close all; clc;
%%
fs = 20000;          
L = 10;           
n = -(L-1)/2 : (L-1)/2;   
n1 = 0:L-1;

f_low = 275;         
f_high = 3250;       

w_low = 2*pi*f_low/fs;
w_high = 2*pi*f_high/fs;

h_ideal = (sin(w_high*n) - sin(w_low*n)) ./ (pi*n);

hamming_window = hamming(L);
blackman_window = blackman(L);
rect = ones(L,1);
triangle = triang(L);

figure
plot([hamming_window, blackman_window, rect, triangle])
legend('Hamming', 'Blackman', 'Rectangular', 'Triangular')
title('Comparison of Window Functions')
h = h_ideal .* rect';

figure;
stem(n1, h, 'filled');
title('Impulse Response');
xlabel('n'); ylabel('h[n]');
grid on;

% Frequency Response
[H, f] = freqz(h, 1, 1024, fs);

% Magnitude Response
figure;
plot(f, abs(H)/max(abs(H)));
title('Magnitude Response (Windowing)');
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
grid on;
xlim([0 fs/2]);  
yline(0.707, 'r--', '-3 dB');
hold off 

% Phase Response
figure;
plot(f, unwrap(angle(H)));
title('Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;

% Pole-Zero Plot
figure;
zplane(roots(h),[]);
title('Pole-Zero Plot');

% Z-transform
syms z
H_z = 0;
for k = 1:L
    H_z = H_z + h(k) * z^(-(k-1));
end

