% Carregar o pacote de processamento de sinais
pkg load signal;

% Carregando os sinais do arquivo sinais.mat
saveVarsMat = load('sinais.mat');
order_butter_filter = 3

x1 = saveVarsMat.x1; % Sinal 1 (amostrado a F1 = 8000 Hz)
x2 = saveVarsMat.x2; % Sinal 2 (amostrado a F2 = 96000 Hz)

% Frequências de amostragem
F1 = 8000; % Frequência de amostragem do sinal x1 (Hz)
F2 = 96000; % Frequência de amostragem do sinal x2 (Hz)

% Análise do sinal x1
T1 = 1/F1; % Período de amostragem do sinal x1
N1 = length(x1); % Número de amostras do sinal x1
t1 = (0:N1-1)*T1; % Vetor de tempo para x1

% Análise do sinal x2
T2 = 1/F2; % Período de amostragem do sinal x2
N2 = length(x2); % Número de amostras do sinal x2
t2 = (0:N2-1)*T2; % Vetor de tempo para x2

%% Figura 1: Plot dos sinais originais e seus espectros
figure;

% Sinal x1 no domínio do tempo
subplot(2,2,1);
plot(t1, x1);
title('Sinal x1 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Sinal x2 no domínio do tempo
subplot(2,2,2);
plot(t2, x2);
title('Sinal x2 no Domínio do Tempo');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Cálculo do espectro de frequências usando a FFT
X1 = fft(x1); % Transformada de Fourier do sinal x1
X1_magnitude = abs(X1); % Magnitude do espectro
f1 = (0:N1-1)*(F1/N1); % Vetor de frequências para x1

X2 = fft(x2); % Transformada de Fourier do sinal x2
X2_magnitude = abs(X2); % Magnitude do espectro
f2 = (0:N2-1)*(F2/N2); % Vetor de frequências para x2

% Plot dos espectros de frequências dos sinais originais
subplot(2,2,3);
plot(f1, X1_magnitude);
title('Espectro de Frequências do Sinal x1');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,2,4);
plot(f2, X2_magnitude);
title('Espectro de Frequências do Sinal x2');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

%% Figura 2: Manipulações no sinal x1 (filtragem, normalização, etc.)
figure;

% Plot do sinal x1 no domínio do tempo (antes do filtro)
subplot(3,2,1);
plot(t1, x1);
title('Sinal x1 no Domínio do Tempo (Antes do Filtro)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Plot do espectro de frequências do sinal x1 (antes do filtro)
subplot(3,2,2);
plot(f1, X1_magnitude);
title('Espectro de Frequências do Sinal x1 (Antes do Filtro)');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Projeto do filtro anti-aliasing (passa-baixa) para x1
fc1 = F1/2; % Frequência de corte do filtro (metade da frequência de amostragem)
[b1, a1] = butter(order_butter_filter, fc1/(F1/2)); % Filtro Butterworth de 3ª ordem

% Aplicação do filtro anti-aliasing ao sinal x1
x1_filtered = filter(b1, a1, x1); % Sinal x1 após o filtro

% Plot do sinal x1 no domínio do tempo (após o filtro)
subplot(3,2,3);
plot(t1, x1_filtered);
title('Sinal x1 no Domínio do Tempo (Após o Filtro)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Cálculo do espectro de frequências usando a FFT (após o filtro)
X1_filtered = fft(x1_filtered); % Transformada de Fourier do sinal x1 filtrado
X1_filtered_magnitude = abs(X1_filtered); % Magnitude do espectro

% Plot do espectro de frequências do sinal x1 (após o filtro)
subplot(3,2,4);
plot(f1, X1_filtered_magnitude);
title('Espectro de Frequências do Sinal x1 (Após o Filtro)');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Verificação do Teorema da Amostragem após o filtro (x1)
f_max_x1_filtered = max(f1(X1_filtered_magnitude > 0.1*max(X1_filtered_magnitude))); % Frequência máxima significativa
if F1 < 2*f_max_x1_filtered
    warning('Atenção: Frequência de amostragem de x1 viola o Teorema de Nyquist! Pode ocorrer aliasing.');
else
    disp('Frequência de amostragem de x1 adequada conforme o Teorema de Nyquist.');
end

% Preparando o sinal x1 filtrado para manipulação
% Exemplo: Normalização do sinal
x1_filtered_normalized = x1_filtered / max(abs(x1_filtered)); % Normalização do sinal

% Verificação da normalização (x1)
if any(abs(x1_filtered_normalized) > 1)
    warning('Atenção: O sinal x1 normalizado contém valores fora do intervalo [-1, 1].');
else
    disp('Sinal x1 normalizado corretamente no intervalo [-1, 1].');
end

% Plot do sinal x1 filtrado e normalizado (pronto para manipulação)
subplot(3,2,5);
plot(t1, x1_filtered_normalized);
title('Sinal x1 Filtrado e Normalizado (Pronto para Manipulação)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Salvando o sinal x1 filtrado e normalizado em um arquivo WAV na pasta raiz
output_filename_x1 = './sinal_filtrado_x1.wav'; % Caminho relativo para a pasta raiz
try
    audiowrite(output_filename_x1, x1_filtered_normalized, F1); % Salva o sinal no formato WAV
    disp(['Sinal x1 filtrado salvo como: ', output_filename_x1]);
    disp(['Caminho completo: ', fullfile(pwd, output_filename_x1)]);
catch ME
    warning('Erro ao salvar o arquivo WAV (x1):');
    disp(ME.message);
end

%% Figura 3: Manipulações no sinal x2 (filtragem, normalização, etc.)
figure;

% Plot do sinal x2 no domínio do tempo (antes do filtro)
subplot(3,2,1);
plot(t2, x2);
title('Sinal x2 no Domínio do Tempo (Antes do Filtro)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Plot do espectro de frequências do sinal x2 (antes do filtro)
subplot(3,2,2);
plot(f2, X2_magnitude);
title('Espectro de Frequências do Sinal x2 (Antes do Filtro)');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Projeto do filtro anti-aliasing (passa-baixa) para x2
fc2 = F2/2; % Frequência de corte do filtro (metade da frequência de amostragem)
[b2, a2] = butter(order_butter_filter, fc2/(F2/2)); % Filtro Butterworth de 4ª ordem

% Aplicação do filtro anti-aliasing ao sinal x2
x2_filtered = filter(b2, a2, x2); % Sinal x2 após o filtro

% Plot do sinal x2 no domínio do tempo (após o filtro)
subplot(3,2,3);
plot(t2, x2_filtered);
title('Sinal x2 no Domínio do Tempo (Após o Filtro)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Cálculo do espectro de frequências usando a FFT (após o filtro)
X2_filtered = fft(x2_filtered); % Transformada de Fourier do sinal x2 filtrado
X2_filtered_magnitude = abs(X2_filtered); % Magnitude do espectro

% Plot do espectro de frequências do sinal x2 (após o filtro)
subplot(3,2,4);
plot(f2, X2_filtered_magnitude);
title('Espectro de Frequências do Sinal x2 (Após o Filtro)');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Verificação do Teorema da Amostragem após o filtro (x2)
f_max_x2_filtered = max(f2(X2_filtered_magnitude > 0.1*max(X2_filtered_magnitude))); % Frequência máxima significativa
if F2 < 2*f_max_x2_filtered
    warning('Atenção: Frequência de amostragem de x2 viola o Teorema de Nyquist! Pode ocorrer aliasing.');
else
    disp('Frequência de amostragem de x2 adequada conforme o Teorema de Nyquist.');
end

% Preparando o sinal x2 filtrado para manipulação
% Exemplo: Normalização do sinal
x2_filtered_normalized = x2_filtered / max(abs(x2_filtered)); % Normalização do sinal

% Verificação da normalização (x2)
if any(abs(x2_filtered_normalized) > 1)
    warning('Atenção: O sinal x2 normalizado contém valores fora do intervalo [-1, 1].');
else
    disp('Sinal x2 normalizado corretamente no intervalo [-1, 1].');
end

% Plot do sinal x2 filtrado e normalizado (pronto para manipulação)
subplot(3,2,5);
plot(t2, x2_filtered_normalized);
title('Sinal x2 Filtrado e Normalizado (Pronto para Manipulação)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Salvando o sinal x2 filtrado e normalizado em um arquivo WAV na pasta raiz
output_filename_x2 = './sinal_filtrado_x2.wav'; % Caminho relativo para a pasta raiz
try
    audiowrite(output_filename_x2, x2_filtered_normalized, F2); % Salva o sinal no formato WAV
    disp(['Sinal x2 filtrado salvo como: ', output_filename_x2]);
    disp(['Caminho completo: ', fullfile(pwd, output_filename_x2)]);
catch ME
    warning('Erro ao salvar o arquivo WAV (x2):');
    disp(ME.message);
end
