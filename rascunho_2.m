% Carregar o pacote de processamento de sinais
pkg load signal;

% Carregando os sinais do arquivo sinais.mat
saveVarsMat = load('sinais.mat');

x1 = saveVarsMat.x1; % Sinal 1 (amostrado a F1 = 8000 Hz)
x2 = saveVarsMat.x2; % Sinal 2 (amostrado a F2 = 96000 Hz)

% Frequências de amostragem
F1 = 8000; % Frequência de amostragem do sinal x1 (Hz)
F2 = 96000; % Frequência de amostragem do sinal x2 (Hz)

% Análise do sinal x1
T1 = 1/F1; % Período de amostragem do sinal x1
N1 = length(x1); % Número de amostras do sinal x1
t1 = (0:N1-1)*T1; % Vetor de tempo para x1

% Plot do sinal x1 no domínio do tempo (antes do filtro)
figure;
subplot(3,2,1);
plot(t1, x1);
title('Sinal x1 no Domínio do Tempo (Antes do Filtro)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Cálculo do espectro de frequências usando a FFT (antes do filtro)
X1 = fft(x1); % Transformada de Fourier do sinal x1
X1_magnitude = abs(X1); % Magnitude do espectro
f1 = (0:N1-1)*(F1/N1); % Vetor de frequências para x1

% Plot do espectro de frequências do sinal x1 (antes do filtro)
subplot(3,2,2);
plot(f1, X1_magnitude);
title('Espectro de Frequências do Sinal x1 (Antes do Filtro)');
xlabel('Frequência (Hz)');
ylabel('Magnitude');
grid on;

% Projeto do filtro anti-aliasing (passa-baixa)
fc = F1/2; % Frequência de corte do filtro (metade da frequência de amostragem)
[b, a] = butter(3, fc/(F1/2)); % Filtro Butterworth de 4ª ordem

% Aplicação do filtro anti-aliasing
x1_filtered = filter(b, a, x1); % Sinal x1 após o filtro

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

% Verificação do Teorema da Amostragem após o filtro
f_max_x1_filtered = max(f1(X1_filtered_magnitude > 0.1*max(X1_filtered_magnitude))); % Frequência máxima significativa
if F1 < 2*f_max_x1_filtered
    warning('Atenção: Frequência de amostragem de x1 viola o Teorema de Nyquist! Pode ocorrer aliasing.');
else
    disp('Frequência de amostragem de x1 adequada conforme o Teorema de Nyquist.');
end

% Preparando o sinal x1 filtrado para manipulação
% Exemplo: Normalização do sinal
x1_filtered_normalized = x1_filtered / max(abs(x1_filtered)); % Normalização do sinal

% Verificação da normalização
if any(abs(x1_filtered_normalized) > 1)
    warning('Atenção: O sinal normalizado contém valores fora do intervalo [-1, 1].');
else
    disp('Sinal normalizado corretamente no intervalo [-1, 1].');
end

% Plot do sinal x1 filtrado e normalizado (pronto para manipulação)
subplot(3,2,5);
plot(t1, x1_filtered_normalized);
title('Sinal x1 Filtrado e Normalizado (Pronto para Manipulação)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Salvando o sinal filtrado e normalizado em um arquivo WAV na pasta raiz
output_filename = './sinal_filtrado.wav'; % Caminho relativo para a pasta raiz
try
    audiowrite(output_filename, x1_filtered_normalized, F1); % Salva o sinal no formato WAV
    disp(['Sinal filtrado salvo como: ', output_filename]);
    disp(['Caminho completo: ', fullfile(pwd, output_filename)]);
catch ME
    warning('Erro ao salvar o arquivo WAV:');
    disp(ME.message);
end
