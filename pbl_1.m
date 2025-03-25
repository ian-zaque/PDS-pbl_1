% Carregar o pacote de processamento de sinais
pkg load signal;

% Carregando os sinais do arquivo sinais.mat
saveVarsMat = load('sinais.mat');

x1 = saveVarsMat.x1; % Sinal 1 (amostrado a F1 = 8000 Hz)
x2 = saveVarsMat.x2; % Sinal 2 (amostrado a F2 = 96000 Hz)
clear saveVarsMat;

% Frequências e Períodos de amostragem
Fs_1 = 8000; % Frequência de amostragem do sinal x1 (Hz)
Ts_1 = 1/(8000); % Periodo de amostragem do sinal x1 (s)
Fs_2 = 96000; % Frequência de amostragem do sinal x2 (Hz)
Ts_2 = 1/(96000); % Periodo de amostragem do sinal x2 (s)

lp_filter_freq_1 = Fs_1/2;
lp_filter_freq_2 = Fs_2/2;
order_butter_filter = 3;

% Sinais x1 e x2 no domínio do tempo
N1 = length(x1); % Número de amostras do sinal x1
xc_1 = (0: N1 - 1) * Ts_1; % Vetor de tempo para x1

N2 = length(x2); % Número de amostras do sinal x2
xc_2 = (0: N2 - 1) * Ts_2; % Vetor de tempo para x2

%Figura 1: Plot dos sinais originais no domínio do tempo
figure;

% Sinal x1 no domínio do tempo
subplot(2,2,1);
plot(xc_1, x1);
title('Sinal x1 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Sinal x2 no domínio do tempo
subplot(2,2,2);
plot(xc_2, x2);
title('Sinal x2 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Do dominio do tempo para frequencia
% Assim, calcularemos o espectro de frequencias de x1 e x2
X1 = fft(x1); % Transformada rapida de fourier de x1
X1_magnitude = abs(X1); % Magnitude do espectro
F_x1 = (0: N1 - 1) * (Fs_1/N1); % Vetor de frequências para x1

X2 = fft(x2); % Transformada rapida de fourier de x2
X2_magnitude = abs(X2); % Magnitude do espectro
F_x2 = (0: N2 - 1) * (Fs_2/N2); % Vetor de frequências para x2

% Plot dos espectros de frequências dos sinais originais x1 e x2
subplot(2,2,3);
plot(F_x1, X1_magnitude);
title('Espectro de Frequências do Sinal x1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,4);
plot(F_x2, X2_magnitude);
title('Espectro de Frequências do Sinal x2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Normalização de X1 e X2
% Removendo repetição do sinal transformado
X1_normalized = abs(X1) / N1;
X2_normalized = abs(X2) / N2;

F_X1_norm = (0:N1 - 1) * (Fs_1/N1);
F_X2_norm = (0:N2 - 1) * (Fs_2/N2);

figure;
subplot(3,2,1);
plot(F_X1_norm(1:N1/2), X1_normalized(1:N1/2));
title('Espectro Normalizado de Frequências do Sinal X1');
xlabel('Frequencia (hz)');
ylabel('Magnitude');
grid on;

subplot(3,2,2);
plot(F_X2_norm(1:N2/2), X2_normalized(1:N2/2));
title('Espectro Normalizado de Frequências do Sinal X2');
xlabel('Frequencia (hz)');
ylabel('Magnitude');
grid on;

% Filtro passa baixa em x1 e x2
%[b1, a1] = butter(order_butter_filter, lp_filter_freq_1/(Fs_1/2));
%x1_filtered = filter(b1, a1, x1);

%[b2, a2] = butter(order_butter_filter-1, lp_filter_freq_2/(Fs_2/2));
%x2_filtered = filter(b2, a2, x2);

%figure;
%subplot(3,2,1);
%plot(xc_1, x1_filtered);
%title('Sinal x1 filtrado no dominio do tempo');
%xlabel('Tempo (s)');
%ylabel('Magnitude');
%grid on;

%subplot(3,2,2);
%plot(xc_2, x2_filtered);
%title('Sinal x2 filtrado no dominio do tempo');
%xlabel('Tempo (s)');
%ylabel('Magnitude');
%grid on;

% Do dominio do tempo para frequencia
% Assim, calcularemos o espectro de frequencias de x1_filtered e x2_filtered
%X1_filtered = fft(x1_filtered); % Transformada rapida de fourier de x1
%X1_filtered_magnitude = abs(X1_filtered); % Magnitude do espectro
%F_filtered_x1 = (0: N1 - 1) * (Fs_1/N1); % Vetor de frequências para x1_filtered

%subplot(3,2,3);
%plot(F_filtered_x1, X1_filtered_magnitude);
%title('Sinal x1 filtrado no dominio da frequencia');
%xlabel('Frequencia (hz)');
%ylabel('Magnitude');
%grid on;
