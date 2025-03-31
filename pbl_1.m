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

%Variaveis de diretorio
currentDir = pwd;
file_path = fileparts(mfilename('fullpath'));
audioDir = fullfile(file_path, 'audios');
x1_filename = 'x1_resample_48k.wav';
x2_filename = 'x2_resample_48k.wav';

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
fs_resample = 48000; % Frequencia de reamostragem comum para x1 e x2 (Hz)
order_butter_filter = 3;

% Sinais x1 e x2 no domínio do tempo
x1_length = length(x1); % Número de amostras do sinal x1
x2_length = length(x2); % Número de amostras do sinal x2
t_x1 = (0: x1_length - 1) / fs_1; % Vetor de tempo para x1
t_x2 = (0: x2_length - 1) / fs_2; % Vetor de tempo para x2

% Do dominio do tempo para frequencia
% Assim, calcularemos o espectro de frequencias de x1 e x2
% Transformadas de Fourier para x1 e x2
X1 = fft(x1); % Transformada rapida de fourier de x1
X1 = fftshift(X1);  % Centralizando (deslocamento)
X1_magnitude = abs(X1);
f_x1 = linspace(-fs_1/2, fs_1/2, x1_length); % Vetor de frequências para x1

X2 = fft(x2); % Transformada rapida de fourier de x2
X2 = fftshift(X2);  % Centralizando (deslocamento)
X2_magnitude = abs(X2);
f_x2 = linspace(-fs_2/2, fs_2/2, x2_length); % Vetor de frequências para x2

%Figura 1: Plot dos sinais originais no domínio do tempo
figure;

% Sinal x1 no domínio do tempo
subplot(3,2,1);
plot(t_x1, x1);
title('Sinal x1 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Sinal x2 no domínio do tempo
subplot(3,2,2);
plot(t_x2, x2);
title('Sinal x2 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Espectros de frequências do sinal original x1
subplot(3,2,3);
plot(f_x1, X1_magnitude);
title('Espectro de Frequências do Sinal x1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Espectros de frequências do sinal original x2
subplot(3,2,4);
plot(f_x2, X2_magnitude);
title('Espectro de Frequências do Sinal x2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;


% Reamostragem do espectro dos sinais X1 e X2
% X1: de 8khz para 48khz (upsample)    | Aumentar 0s entre amostras em 6
% X2: de 96khz para 48khz (downsample) | Decimar amostras em 2
% Criando e aplicando anti-aliasing em x2 antes da decimação para evitar aliasing
[b, a] = butter(order_butter_filter, (fs_resample/fs_2) * 0.9);
x2_filtered = filter(b, a, x2);

x1_up = resample(x1, 6, 1); % Upsampling de fs de 8 kHz para 48 kHz (6 para 1)
x2_down = resample(x2_filtered, 1, 2); % Downsampling de fs 96 kHz para 48 kHz (1 para 2)

x1_up_length = length(x1_up);
x2_down_length = length(x2_down);
t_x1_up = (0: x1_up_length - 1) / fs_resample; % Vetor de tempo para x1_up
t_x2_down = (0: x2_down_length - 1) / fs_resample; % Vetor de tempo para x2_down

X1_up = fft(x1_up); % Transformada rapida de fourier de x1_up
X1_up = fftshift(X1_up);  % Centralizando (deslocamento)
X2_down = fft(x2_down); % Transformada rapida de fourier de x2_down
X2_down = fftshift(X2_down);  % Centralizando (deslocamento)

% Vetores de frequência reamostrados
X1_up_length = length(X1_up);
f_X1_up = linspace(-fs_resample/2, fs_resample/2, X1_up_length);
X2_down_length = length(X2_down);
f_X2_down = linspace(-fs_resample/2, fs_resample/2, X2_down_length);

%Figura 2: Plot dos sinais reamostrados no domínio do tempo
figure;

% Sinal reamostrado x1 no domínio do tempo
subplot(3,2,1);
plot(t_x1_up, x1_up);
title('Sinal x1_up no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Sinal reamostrado x2 no domínio do tempo
subplot(3,2,2);
plot(t_x2_down, x2_down);
title('Sinal x2_down no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Espectros de frequências do sinal reamostrado x1
subplot(3,2,3);
plot(f_X1_up, abs(X1_up));
title('Espectro de Frequências do Sinal x1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Espectros de frequências do sinal reamostrado x2
subplot(3,2,4);
plot(f_X2_down, abs(X2_down));
title('Espectro de Frequências do Sinal x2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

if ~exist(audioDir, 'dir')
    mkdir(audioDir);
end

output_filename_X1 = fullfile(audioDir, x1_filename);
output_filename_X2 = fullfile(audioDir, x2_filename);

audiowrite(output_filename_X1, x1_up, fs_resample);
audiowrite(output_filename_X2, x2_down, fs_resample);

