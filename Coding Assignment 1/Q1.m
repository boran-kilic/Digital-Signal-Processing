close all
clear 
clc

rng(2,"twister")
N = 32;
real_part = randn(1,N);
imag_part = randn(1,N);
x = real_part + 1i*imag_part;
steps_a_h(x);

rng(2,"twister")
N = 256;
real_part = randn(1,N);
imag_part = randn(1,N);
x = real_part + 1i*imag_part;
steps_a_h(x);

rng(2,"twister")
N = 2^12;
real_part = randn(1,N);
imag_part = randn(1,N);
x = real_part + 1i*imag_part;
steps_a_h(x);

function steps_a_h(x)
N = length(x);
fprintf('N = %d \n', N);
n = 0:N-1; 

figure;
subplot(2,1,1);
stem(n, real(x), 'filled');
title('Real Part of x[n]');
xlabel('n');
ylabel('Re\{x[n]\}');
xlim([0 N-1]);
grid on;

subplot(2,1,2);
stem(n, imag(x), 'filled');
title('Imaginary Part of x[n]');
xlabel('n');
ylabel('Im\{x[n]\}');
xlim([0 N-1]);
grid on;

t1 = tic;
Xsum = DFT_summation(x);
elapsed_time = toc(t1);
fprintf('Elapsed time for DFT_summation is %.6f seconds.\n', elapsed_time);
plotDFT(Xsum, 'Summation Formula')

t2 = tic;
X_mat = DFT_matrix(x);
elapsed_time = toc(t2);
fprintf('Elapsed time for DFT_matrix is %.6f seconds.\n', elapsed_time);
plotDFT(X_mat, 'DFT matrix')

t3 = tic;
X_DIT = FFT_DIT(x);
elapsed_time = toc(t3);
fprintf('Elapsed time for FFT_DIT is %.6f seconds.\n', elapsed_time);
plotDFT(X_DIT, 'FFT-DIT')

t4 = tic;
X_DIF = bitrevorder(FFT_DIF(x));
elapsed_time = toc(t4);
fprintf('Elapsed time for FFT_DIF is %.6f seconds.\n', elapsed_time);
plotDFT(X_DIF, 'FFT-DIF')

t5 = tic;
X_fft = fft(x);
elapsed_time = toc(t5);
fprintf('Elapsed time for fft is %.6f seconds.\n', elapsed_time);
plotDFT(X_fft, 'FFT')
fprintf('\n');

disp(['Summation Formula:   ', num2str(norm(X_fft - Xsum))]);
disp(['DFT matrix:          ', num2str(norm(X_fft - X_mat))]);
disp(['FFT-DIT:             ', num2str(norm(X_fft - X_DIT))]);
disp(['FFT-DIF:             ', num2str(norm(X_fft - X_DIF))]);
disp(['fft:                 ', num2str(norm(X_fft - X_fft))]);
fprintf('\n');
end

function plotDFT(X, str)
N = length(X);
n = 0:N-1;

figure;
subplot(2,1,1);
stem(n, abs(X), 'filled');
title(['Magnitude of X[k] using ', str]);
xlabel('k');
ylabel('|X[k]|');
xlim([0 N-1]);
grid on;

subplot(2,1,2);
stem(n, angle(X), 'filled');
title(['Phase of X[k] using ', str]);
xlabel('k');
ylabel('∠X[k] (radians)');
xlim([0 N-1]);
grid on;
end

function X = FFT_DIF(x)
    N = length(x);
    if N == 1
        X = x;
    else
        X = x;
        % Butterfly stage
        for k = 1:N/2
            temp = X(k);
            X(k) = temp + X(k + N/2);
            X(k + N/2) = (temp - X(k + N/2)) * exp(-1i * 2 * pi * (k - 1) / N);
        end
        % Recursive stage
        X(1:N/2) = FFT_DIF(X(1:N/2));
        X(N/2+1:N) = FFT_DIF(X(N/2+1:N)); 
    end
end


function X = FFT_DIT(x)
    N = length(x);
    if N == 1
        X = x;
    else
        % Divide
        x_even = x(1:2:end);
        x_odd = x(2:2:end);
        
        % Conquer
        X_even = FFT_DIT(x_even);
        X_odd = FFT_DIT(x_odd);
        
        % Combine
        WN = exp(-1i * 2 * pi / N);
        W = 1;
        X = zeros(1, N);
        for k = 1:N/2
            X(k) = X_even(k) + W * X_odd(k);
            X(k + N/2) = X_even(k) - W * X_odd(k);
            W = W * WN;
        end
    end
end


function X = DFT_matrix(x)
N = length(x);
WN = exp(-1i * 2 * pi / N); 
DFT_mat = zeros(N, N); 
for k = 0:N-1
    for n = 0:N-1
        DFT_mat(k+1, n+1) = WN^(k * n);
    end
end
X = (DFT_mat * x.').';
end


function X = DFT_summation(x)
N = length(x);
X = zeros(1, N); 
for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1) * exp(-1i * 2 * pi * k * n / N);
    end
end
end
