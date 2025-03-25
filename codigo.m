% Carregar os dados dos arquivos (assumindo que estão no formato correto)

saveVarsMat = load('sinais.mat');

x1_n = saveVarsMat.x1; % Sinal 1 (amostrado a F1 = 8000 Hz)
x2_n = saveVarsMat.x2;

%x1_n = load('x1_n.txt');  % Carregar os dados de x1[n]
%x2_n = load('x2_n.txt');  % Carregar os dados de x2[n]

% Definir as taxas de amostragem (em Hz)
fs_x1 = 8000;   % Taxa de amostragem original para x1[n] em Hz
fs_x2 = 96000;  % Taxa de amostragem original para x2[n] em Hz

%% Processamento de x2 (Subamostragem)
% Calcular a Transformada Discreta de Fourier (DFT) de x2[n]
X2_k = fft(x2_n);
f2 = (0:length(x2_n)-1)*(fs_x2/length(x2_n));  % Frequências para x2[n] em Hz

% Implementação manual do filtro passa-baixa com frequência de corte de 4 kHz (para D=6)
cutoff_x2 = 4000;  % Frequência de corte em Hz
filter_x2 = zeros(size(f2));  % Inicializa o filtro com zeros

% Criar manualmente o filtro passa-baixa ideal
for i = 1:length(f2)
    if f2(i) <= cutoff_x2 || f2(i) >= (fs_x2 - cutoff_x2)
        filter_x2(i) = 1;  % Permite as frequências abaixo do corte ou as frequências espelhadas
    else
        filter_x2(i) = 0;  % Bloqueia as frequências acima do corte
    end
end

% Aplicar o filtro no domínio da frequência
X2_k_filtered = X2_k .* filter_x2';
x2_n_filtered = ifft(X2_k_filtered);  % Sinal filtrado no domínio do tempo

% Fator de decimação para reduzir de 96 kHz para 16 kHz
D = 6;  % Fator de decimação (96 kHz / 16 kHz = 6)

% Substituição da função downsample por implementação manual
% x2_n_decimated = downsample(real(x2_n_filtered), D);
x2_n_decimated = manual_downsample(real(x2_n_filtered), D);  % Subamostrando o sinal
x2_n_decimated = double(x2_n_decimated);  % Garantir tipo double

% Calcular o espectro do sinal subamostrado
fs_x2_new = fs_x2 / D;  % Nova taxa de amostragem: 16 kHz
X2_k_decimated = fft(x2_n_decimated);
f2_decimated = (0:length(x2_n_decimated)-1)*(fs_x2_new/length(x2_n_decimated));

%% Processamento de x1 (Superamostragem)
% Calcular a Transformada Discreta de Fourier (DFT) de x1[n]
X1_k = fft(x1_n);
f1 = (0:length(x1_n)-1)*(fs_x1/length(x1_n));  % Frequências para x1[n] em Hz

% Fator de interpolação para aumentar de 8 kHz para 16 kHz
L = 2;  % Fator de interpolação (16 kHz / 8 kHz = 2)

% Substituição da função upsample por implementação manual
% x1_n_upsampled = upsample(x1_n, L);
x1_n_upsampled = manual_upsample(x1_n, L);  % Inserir zeros (upsampling)

% Criar o filtro passa-baixa com frequência de corte de 4 kHz
fs_x1_new = fs_x1 * L;  % Nova taxa de amostragem: 16 kHz
f1_upsampled = (0:length(x1_n_upsampled)-1)*(fs_x1_new/length(x1_n_upsampled));
cutoff_x1 = 4000;  % Frequência de corte em Hz (frequência máxima original)

% Implementação manual do filtro passa-baixa ideal
filter_x1 = zeros(size(f1_upsampled));  % Inicializa o filtro com zeros

% Criar manualmente o filtro passa-baixa ideal
for i = 1:length(f1_upsampled)
    if f1_upsampled(i) <= cutoff_x1 || f1_upsampled(i) >= (fs_x1_new - cutoff_x1)
        filter_x1(i) = 1;  % Permite as frequências abaixo do corte ou as frequências espelhadas
    else
        filter_x1(i) = 0;  % Bloqueia as frequências acima do corte
    end
end

% Aplicar o filtro no domínio da frequência
X1_k_upsampled = fft(x1_n_upsampled);
X1_k_filtered = X1_k_upsampled .* filter_x1';
x1_n_interpolated = ifft(X1_k_filtered);  % Sinal interpolado no domínio do tempo
x1_n_interpolated = real(x1_n_interpolated);  % Garantir que seja real
x1_n_interpolated = double(x1_n_interpolated);  % Garantir tipo double

%% Soma dos sinais após processamento multitaxa
% Ajustar o comprimento dos sinais (se necessário)
min_length = min(length(x1_n_interpolated), length(x2_n_decimated));
x1_adjusted = x1_n_interpolated(1:min_length);
x2_adjusted = x2_n_decimated(1:min_length);

% Normalizar os sinais para evitar clipping
x1_normalized = x1_adjusted / max(abs(x1_adjusted));
x2_normalized = x2_adjusted / max(abs(x2_adjusted));

% Soma dos sinais normalizados
x_soma = x1_normalized + x2_normalized;

% Normalizar novamente a soma para evitar clipping
x_soma_normalized = x_soma / max(abs(x_soma));

%% Reproduzir o áudio somado
soundsc(x_soma_normalized, fs_x1_new);  % Reproduz o sinal somado a 16 kHz

%% Plotar os sinais no domínio do tempo
figure('Name', 'Sinais no Domínio do Tempo');
% x1 interpolado
subplot(3,1,1);
t1 = (0:length(x1_normalized)-1)/fs_x1_new;
plot(t1, x1_normalized);
title('x1[n] Superamostrado (16 kHz)');
xlabel('Tempo (s)');
ylabel('Amplitude');

% x2 decimado
subplot(3,1,2);
t2 = (0:length(x2_normalized)-1)/fs_x2_new;
plot(t2, x2_normalized);
title('x2[n] Subamostrado (16 kHz)');
xlabel('Tempo (s)');
ylabel('Amplitude');

% Soma dos sinais
subplot(3,1,3);
t_soma = (0:length(x_soma_normalized)-1)/fs_x1_new;
plot(t_soma, x_soma_normalized);
title('Soma dos Sinais (16 kHz)');
xlabel('Tempo (s)');
ylabel('Amplitude');

sgtitle('Sinais Processados e Somados');

%% Plotar domínios de frequência do sinal somado com replicação
figure('Name', 'Espectro de Frequência da Soma com Replicação');
X_soma = fft(x_soma_normalized);
N = length(x_soma_normalized);
% Reorganizar o espectro para centralizar em zero
X_soma_shifted = fftshift(X_soma);
% Criar eixo de frequência centralizado em 0
f_soma_shifted = (-N/2:N/2-1)*(fs_x1_new/N);

% Parâmetros para replicação do espectro
num_replicas = 3; % Número de réplicas para cada lado
fs = fs_x1_new; % Frequência de amostragem
max_freq_plot = (num_replicas+0.5)*fs; % Limite máximo para a frequência plotada

% Criar eixo de frequência estendido
num_points = (2*num_replicas+1)*N;
f_extended = linspace(-max_freq_plot, max_freq_plot, num_points);

% Criar espectro replicado
X_extended = zeros(size(f_extended));

% Extrair magnitudes do espectro centralizado
X_mag = abs(X_soma_shifted);

% Preencher a parte central (espectro original)
mid_start = floor(num_points/2) - floor(N/2) + 1;
mid_end = mid_start + N - 1;
X_extended(mid_start:mid_end) = X_mag;

% Replicar o espectro para esquerda e direita
for i = 1:num_replicas
    % Réplica à esquerda
    left_start = mid_start - i*N;
    left_end = left_start + N - 1;
    X_extended(left_start:left_end) = X_mag;

    % Réplica à direita
    right_start = mid_end + 1 + (i-1)*N;
    right_end = right_start + N - 1;
    X_extended(right_start:right_end) = X_mag;
end

% Plotar o espectro replicado
plot(f_extended, X_extended);
title('Domínio da Frequência da Soma dos Sinais com Replicação');
xlabel('Frequência (Hz)');
ylabel('|X(f)|');
xlim([-max_freq_plot, max_freq_plot]); % Limitar a visualização
grid on;

%% Plotar domínios de frequência de x2 com replicação
figure('Name', 'Espectro de Frequência de x2[n] com Replicação');

% Original
subplot(2,1,1);
% Reorganizar o espectro para centralizar em zero
X2_k_shifted = fftshift(X2_k);
N2 = length(X2_k);
% Criar eixo de frequência centralizado em 0
f2_shifted = (-N2/2:N2/2-1)*(fs_x2/N2);

% Parâmetros para replicação
num_replicas_x2 = 3;
fs_x2_full = fs_x2; % Frequência de amostragem
max_freq_plot_x2 = (num_replicas_x2+0.5)*fs_x2_full; % Limite máximo

% Criar eixo de frequência estendido
num_points_x2 = (2*num_replicas_x2+1)*N2;
f2_extended = linspace(-max_freq_plot_x2, max_freq_plot_x2, num_points_x2);

% Criar espectro replicado
X2_extended = zeros(size(f2_extended));

% Extrair magnitudes do espectro centralizado
X2_mag = abs(X2_k_shifted);

% Preencher a parte central (espectro original)
mid_start_x2 = floor(num_points_x2/2) - floor(N2/2) + 1;
mid_end_x2 = mid_start_x2 + N2 - 1;
X2_extended(mid_start_x2:mid_end_x2) = X2_mag;

% Replicar o espectro para esquerda e direita
for i = 1:num_replicas_x2
    % Réplica à esquerda
    left_start = mid_start_x2 - i*N2;
    left_end = left_start + N2 - 1;
    X2_extended(left_start:left_end) = X2_mag;

    % Réplica à direita
    right_start = mid_end_x2 + 1 + (i-1)*N2;
    right_end = right_start + N2 - 1;
    X2_extended(right_start:right_end) = X2_mag;
end

% Plotar
plot(f2_extended, X2_extended);
title('Domínio da Frequência de x2[n] Original com Replicação');
xlabel('Frequência (Hz)');
ylabel('|X2(f)|');
xlim([-max_freq_plot_x2, max_freq_plot_x2]);
grid on;

% Subamostrado
subplot(2,1,2);
% Reorganizar o espectro para centralizar em zero
X2d_shifted = fftshift(X2_k_decimated);
N2d = length(X2_k_decimated);
% Criar eixo de frequência centralizado em 0
f2d_shifted = (-N2d/2:N2d/2-1)*(fs_x2_new/N2d);

% Parâmetros para replicação
num_replicas_x2d = 3;
fs_x2d_full = fs_x2_new; % Frequência de amostragem
max_freq_plot_x2d = (num_replicas_x2d+0.5)*fs_x2d_full; % Limite máximo

% Criar eixo de frequência estendido
num_points_x2d = (2*num_replicas_x2d+1)*N2d;
f2d_extended = linspace(-max_freq_plot_x2d, max_freq_plot_x2d, num_points_x2d);

% Criar espectro replicado
X2d_extended = zeros(size(f2d_extended));

% Extrair magnitudes do espectro centralizado
X2d_mag = abs(X2d_shifted);

% Preencher a parte central (espectro original)
mid_start_x2d = floor(num_points_x2d/2) - floor(N2d/2) + 1;
mid_end_x2d = mid_start_x2d + N2d - 1;
X2d_extended(mid_start_x2d:mid_end_x2d) = X2d_mag;

% Replicar o espectro para esquerda e direita
for i = 1:num_replicas_x2d
    % Réplica à esquerda
    left_start = mid_start_x2d - i*N2d;
    left_end = left_start + N2d - 1;
    X2d_extended(left_start:left_end) = X2d_mag;

    % Réplica à direita
    right_start = mid_end_x2d + 1 + (i-1)*N2d;
    right_end = right_start + N2d - 1;
    X2d_extended(right_start:right_end) = X2d_mag;
end

% Plotar
plot(f2d_extended, X2d_extended);
title('Domínio da Frequência de x2[n] Subamostrado (16 kHz) com Replicação');
xlabel('Frequência (Hz)');
ylabel('|X2(f) Subamostrado|');
xlim([-max_freq_plot_x2d, max_freq_plot_x2d]);
grid on;

sgtitle('Espectro de Frequência de x2[n] com Replicação');

%% Plotar domínios de frequência de x1 com replicação
figure('Name', 'Espectro de Frequência de x1[n] com Replicação');

% Original
subplot(2,1,1);
% Reorganizar o espectro para centralizar em zero
X1_k_shifted = fftshift(X1_k);
N1 = length(X1_k);
% Criar eixo de frequência centralizado em 0
f1_shifted = (-N1/2:N1/2-1)*(fs_x1/N1);

% Parâmetros para replicação
num_replicas_x1 = 3;
fs_x1_full = fs_x1; % Frequência de amostragem
max_freq_plot_x1 = (num_replicas_x1+0.5)*fs_x1_full; % Limite máximo

% Criar eixo de frequência estendido
num_points_x1 = (2*num_replicas_x1+1)*N1;
f1_extended = linspace(-max_freq_plot_x1, max_freq_plot_x1, num_points_x1);

% Criar espectro replicado
X1_extended = zeros(size(f1_extended));

% Extrair magnitudes do espectro centralizado
X1_mag = abs(X1_k_shifted);

% Preencher a parte central (espectro original)
mid_start_x1 = floor(num_points_x1/2) - floor(N1/2) + 1;
mid_end_x1 = mid_start_x1 + N1 - 1;
X1_extended(mid_start_x1:mid_end_x1) = X1_mag;

% Replicar o espectro para esquerda e direita
for i = 1:num_replicas_x1
    % Réplica à esquerda
    left_start = mid_start_x1 - i*N1;
    left_end = left_start + N1 - 1;
    X1_extended(left_start:left_end) = X1_mag;

    % Réplica à direita
    right_start = mid_end_x1 + 1 + (i-1)*N1;
    right_end = right_start + N1 - 1;
    X1_extended(right_start:right_end) = X1_mag;
end

% Plotar
plot(f1_extended, X1_extended);
title('Domínio da Frequência de x1[n] Original com Replicação');
xlabel('Frequência (Hz)');
ylabel('|X1(f)|');
xlim([-max_freq_plot_x1, max_freq_plot_x1]);
grid on;

% Superamostrado
subplot(2,1,2);
% Reorganizar o espectro para centralizar em zero
X1_k_interpolated = fft(x1_n_interpolated);  % Espectro do sinal interpolado
X1i_shifted = fftshift(X1_k_interpolated);
N1i = length(X1_k_interpolated);
% Criar eixo de frequência centralizado em 0
f1i_shifted = (-N1i/2:N1i/2-1)*(fs_x1_new/N1i);

% Parâmetros para replicação
num_replicas_x1i = 3;
fs_x1i_full = fs_x1_new; % Frequência de amostragem
max_freq_plot_x1i = (num_replicas_x1i+0.5)*fs_x1i_full; % Limite máximo

% Criar eixo de frequência estendido
num_points_x1i = (2*num_replicas_x1i+1)*N1i;
f1i_extended = linspace(-max_freq_plot_x1i, max_freq_plot_x1i, num_points_x1i);

% Criar espectro replicado
X1i_extended = zeros(size(f1i_extended));

% Extrair magnitudes do espectro centralizado
X1i_mag = abs(X1i_shifted);

% Preencher a parte central (espectro original)
mid_start_x1i = floor(num_points_x1i/2) - floor(N1i/2) + 1;
mid_end_x1i = mid_start_x1i + N1i - 1;
X1i_extended(mid_start_x1i:mid_end_x1i) = X1i_mag;

% Replicar o espectro para esquerda e direita
for i = 1:num_replicas_x1i
    % Réplica à esquerda
    left_start = mid_start_x1i - i*N1i;
    left_end = left_start + N1i - 1;
    X1i_extended(left_start:left_end) = X1i_mag;

    % Réplica à direita
    right_start = mid_end_x1i + 1 + (i-1)*N1i;
    right_end = right_start + N1i - 1;
    X1i_extended(right_start:right_end) = X1i_mag;
end

% Plotar
plot(f1i_extended, X1i_extended);
title('Domínio da Frequência de x1[n] Superamostrado (16 kHz) com Replicação');
xlabel('Frequência (Hz)');
ylabel('|X1(f) Interpolado|');
xlim([-max_freq_plot_x1i, max_freq_plot_x1i]);
grid on;

sgtitle('Espectro de Frequência de x1[n] com Replicação');

%% Função manual de downsample
function output = manual_downsample(signal, factor)
    % Implementação manual da função downsample
    % Mantém apenas uma amostra a cada 'factor' amostras
    output = signal(1:factor:end);
end

%% Função manual de upsample
function output = manual_upsample(signal, factor)
    % Implementação manual da função upsample
    % Insere (factor-1) zeros entre cada amostra do sinal original
    signalLength = length(signal);
    output = zeros(signalLength * factor, 1);
    output(1:factor:end) = signal;
end

% Salvar o áudio somado (opcional)
% audiowrite('soma_x1_x2.wav', x_soma_normalized, fs_x1_new);
