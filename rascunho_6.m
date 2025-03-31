% Carregar o pacote de processamento de sinais
pkg load signal;

% Carregando os sinais do arquivo sinais.mat
saveVarsMat = load('sinais.mat');

x1 = saveVarsMat.x1; % Sinal 1 (amostrado a 8000 Hz)
x2 = saveVarsMat.x2; % Sinal 2 (amostrado a 96000 Hz)
clear saveVarsMat;

% Frequências e períodos de amostragem
fs_1 = 8000;
fs_2 = 96000;
fs_resample = 48000;
order_butter_filter = 3;

% Vetores de tempo
x1_length = length(x1);
t_x1 = (0:x1_length - 1) / fs_1;

x2_length = length(x2);
t_x2 = (0:x2_length - 1) / fs_2;

% FFT dos sinais originais
X1 = fftshift(fft(x1));
X2 = fftshift(fft(x2));

% Vetores de frequência corrigidos
f_x1 = linspace(-fs_1/2, fs_1/2, x1_length);
f_x2 = linspace(-fs_2/2, fs_2/2, x2_length);

% Aplicando filtro antes da decimação para evitar aliasing
[b, a] = butter(order_butter_filter, (fs_resample/fs_2) * 0.9);
x2_filtered = filter(b, a, x2);

% Reamostragem correta
x1_up = resample(x1, 6, 1); % Upsampling de 8 kHz para 48 kHz
x2_down = resample(x2_filtered, 1, 2); % Downsampling de 96 kHz para 48 kHz

% FFT dos sinais reamostrados
X1_up = fftshift(fft(x1_up));
X2_down = fftshift(fft(x2_down));

% Vetores de frequência após reamostragem
X1_up_length = length(X1_up);
f_X1_up = linspace(-fs_resample/2, fs_resample/2, X1_up_length);

X2_down_length = length(X2_down);
f_X2_down = linspace(-fs_resample/2, fs_resample/2, X2_down_length);

% Salvar os sinais reamostrados
output_filename_X1 = fullfile('C:\Users\xaven\OneDrive\Documentos\UEFS\Semestres\2025.1\MI - PDS\PDS-pbl_1\audios', 'x1_freq_48k.wav');
output_filename_X2 = fullfile('C:\Users\xaven\OneDrive\Documentos\UEFS\Semestres\2025.1\MI - PDS\PDS-pbl_1\audios', 'x2_freq_48k.wav');
audiowrite(output_filename_X1, x1_up, fs_resample);
audiowrite(output_filename_X2, x2_down, fs_resample);

% Plotagem dos sinais e espectros
figure;
subplot(3,2,1); plot(t_x1, x1); title('Sinal x1 no Tempo'); xlabel('Tempo (s)'); ylabel('Amplitude'); grid on;
subplot(3,2,2); plot(t_x2, x2); title('Sinal x2 no Tempo'); xlabel('Tempo (s)'); ylabel('Amplitude'); grid on;
subplot(3,2,3); plot(f_x1, abs(X1)); title('Espectro de x1'); xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on;
subplot(3,2,4); plot(f_x2, abs(X2)); title('Espectro de x2'); xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on;
subplot(3,2,5); plot(f_X1_up, abs(X1_up)); title('Espectro de x1 Reamostrado'); xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on;
subplot(3,2,6); plot(f_X2_down, abs(X2_down)); title('Espectro de x2 Reamostrado'); xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on;

