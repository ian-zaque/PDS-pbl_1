% SUMÁRIO DAS VARIÁVEIS:
% - saveVarsMat: VARIÁVEL TEMPORÁRIA PARA ARMAZENAR DADOS DO ARQUIVO .MAT
% - x1, x2: SINAIS DE ENTRADA AMOSTRADOS EM DIFERENTES FREQUÊNCIAS (8000 HZ E 96000 HZ, RESPECTIVAMENTE)
% - fs_1, fs_2: FREQUÊNCIAS DE AMOSTRAGEM DOS SINAIS x1 E x2
% - ts_1, ts_2: PERÍODOS DE AMOSTRAGEM CORRESPONDENTES
% - lp_filter_freq_1, lp_filter_freq_2: FREQUÊNCIAS DE CORTE DOS FILTROS PASSA-BAIXA
% - order_butter_filter: ORDEM DO FILTRO BUTTERWORTH
% - x1_length, x2_length: NÚMERO DE AMOSTRAS NOS SINAIS x1 E x2
% - t_x1, t_x2: VETORES DE TEMPO PARA OS SINAIS x1 E x2
% - X1, X2: TRANSFORMADAS DE FOURIER DOS SINAIS x1 E x2
% - X1_magnitude, X2_magnitude: MAGNITUDE DOS ESPECTROS DE FREQUÊNCIA DOS SINAIS
% - f_x1, f_x2: VETORES DE FREQUÊNCIA ASSOCIADOS AOS ESPECTROS DE x1 E x2
% - X1_up, X2_down: TRANSFORMADAS DE FOURIER DOS SINAIS REAMOSTRADOS
% - X1_up_length, X2_down_length: NÚMERO DE AMOSTRAS NOS SINAIS REAMOSTRADOS
% - X1_up_magnitude, X2_down_magnitude: MAGNITUDE DOS ESPECTROS DE FREQUÊNCIA DOS SINAIS REAMOSTRADOS
% - f_X1_up, f_X2_down: VETORES DE FREQUÊNCIA ASSOCIADOS AOS ESPECTROS REAMOSTRADOS


% Carregar o pacote de processamento de sinais
pkg load signal;

% Carregando os sinais do arquivo sinais.mat
saveVarsMat = load('sinais.mat');

x1 = saveVarsMat.x1; % Sinal 1 (amostrado a F1 = 8000 Hz)
x2 = saveVarsMat.x2; % Sinal 2 (amostrado a F2 = 96000 Hz)
clear saveVarsMat;

% Frequências e Períodos de amostragem
fs_1 = 8000; % Frequência de amostragem do sinal x1 (Hz)
fs_2 = 96000; % Frequência de amostragem do sinal x2 (Hz)
ts_1 = 1/(8000); % Periodo de amostragem do sinal x1 (s)
ts_2 = 1/(96000); % Periodo de amostragem do sinal x2 (s)

lp_filter_freq_1 = fs_1/2;
lp_filter_freq_2 = fs_2/2;
order_butter_filter = 3;

% Sinais x1 e x2 no domínio do tempo
x1_length = length(x1); % Número de amostras do sinal x1
t_x1 = (0: x1_length - 1) * ts_1; % Vetor de tempo para x1

x2_length = length(x2); % Número de amostras do sinal x2
t_x2 = (0: x2_length - 1) * ts_2; % Vetor de tempo para x2

%Figura 1: Plot dos sinais originais no domínio do tempo
figure;

% Sinal x1 no domínio do tempo
subplot(2,2,1);
plot(t_x1, x1);
title('Sinal x1 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Sinal x2 no domínio do tempo
subplot(2,2,2);
plot(t_x2, x2);
title('Sinal x2 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Do dominio do tempo para frequencia
% Assim, calcularemos o espectro de frequencias de x1 e x2
X1 = fft(x1); % Transformada rapida de fourier de x1
X1_magnitude = abs(X1); % Magnitude do espectro
f_x1 = (0: x1_length - 1) * (fs_1/x1_length); % Vetor de frequências para x1

X2 = fft(x2); % Transformada rapida de fourier de x2
X2_magnitude = abs(X2); % Magnitude do espectro
f_x2 = (0: x2_length - 1) * (fs_2/x2_length); % Vetor de frequências para x2

% Plot dos espectros de frequências dos sinais originais x1 e x2
subplot(2,2,3);
plot(f_x1, X1_magnitude);
title('Espectro de Frequências do Sinal x1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,4);
plot(f_x2, X2_magnitude);
title('Espectro de Frequências do Sinal x2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Reamostragem do espectro dos sinais X1 e X2
% X1: de 8khz para 16khz (upsample)    | Interpolar em 2
% X2: de 96khz para 16khz (downsample) | Decimar em 6
X1_up = upsample(X1, 2);
X1_up_length = length(X1_up);
X1_up_magnitude = abs(X1_up);
f_X1_up = (0: X1_up_length - 1) * ( (fs_1*2) /X1_up_length); % Vetor de frequências para X1_up

X2_down = downsample(X2, 6);
X2_down_length = length(X2_down);
X2_down_magnitude = abs(X2_down);
f_X2_down = (0: X2_down_length - 1) * ( (fs_2/6) /X2_down_length); % Vetor de frequências para X2_down

% APLICANDO A TRANSFORMADA INVERSA DE FOURIER PARA TER O SINAL NO DOMINIO DO TEMPO
X1_t_up = ifft(X1_up);
X1_t_up_length = length(X1_t_up);
t_X1_up = (0: X1_t_up_length - 1) / (fs_1*2);

X2_t_down = ifft(X2_down);
X2_t_down_length = length(X2_t_down);
t_X2_down = (0: X2_t_down_length - 1) / (fs_2/6);

figure;
subplot(3,2,1);
plot(f_X1_up, X1_up_magnitude);
title('Espectro de Frequências do Sinal x1 Reamostrado');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

subplot(3,2,2);
plot(f_X2_down, X2_down_magnitude);
title('Espectro de Frequências do Sinal x2 Reamostrado');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

subplot(3,2,3);
plot(t_X1_up, X1_t_up);
title('Sinal X1 reamostrado no dominio do tempo');
xlabel('Tempo (s)');
ylabel('Magnitude');
grid on;

subplot(3,2,4);
plot(t_X2_down, X2_t_down);
title('Sinal X2 reamostrado no dominio do tempo');
xlabel('Tempo (s)');
ylabel('Magnitude');
grid on;

% Normalização de X1 e X2
% Removendo repetição do sinal transformado
%X1_normalized = abs(X1) / x1_length;
%X2_normalized = abs(X2) / x2_length;

%f_x1_norm = (0:x1_length - 1) * (fs_1/x1_length);
%f_x2_norm = (0:x2_length - 1) * (fs_2/x2_length);

%figure;
%subplot(4,2,1);
%plot(f_x1_norm(1:x1_length/2), X1_normalized(1:x1_length/2));
%title('Espectro Normalizado de Frequências do Sinal X1');
%xlabel('Frequencia (hz)');
%ylabel('Magnitude');
%grid on;

%subplot(4,2,2);
%plot(f_x2_norm(1:x2_length/2), X2_normalized(1:x2_length/2));
%title('Espectro Normalizado de Frequências do Sinal X2');
%xlabel('Frequencia (hz)');
%ylabel('Magnitude');
%grid on;

