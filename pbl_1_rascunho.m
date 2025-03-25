% Carregar os sinais do arquivo 'sinais.mat'
saveVarsMat = load('sinais.mat');

x1 = saveVarsMat.x1; % Sinal 1
x2 = saveVarsMat.x2; % Sinal 2

% Definir as taxas de amostragem
F1_s = 8000;  % Taxa de amostragem para o sinal 1 em Hz
F2_s = 96000; % Taxa de amostragem para o sinal 2 em Hz

% Número de amostras
N_1 = length(x1);
N_2 = length(x2);

% Aplicar a Transformada Rápida de Fourier (FFT)
X1_f = fft(x1);
X2_f = fft(x2);

% Normalizar a magnitude da FFT
X1_f = fftshift(abs(X1_f) / N_1);  % Normalização e centralização para o sinal 1
X2_f = fftshift(abs(X2_f) / N_2);  % Normalização e centralização para o sinal 2

% Criar vetor de frequências
f_1 = linspace(-F1_s/2, F1_s/2, N_1);  % Frequências para o sinal 1
f_2 = linspace(-F2_s/2, F2_s/2, N_2);  % Frequências para o sinal 2

% Plotar o espectro de frequência para o sinal 1
figure;
plot(f_1, X1_f);
title('Espectro de Frequência do Sinal 1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Plotar o espectro de frequência para o sinal 2
figure;
plot(f_2, X2_f);
title('Espectro de Frequência do Sinal 2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

