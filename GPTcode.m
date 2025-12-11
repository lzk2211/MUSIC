% 参数设置
M = 4; % 天线个数
d = 0.5; % 天线间距（以波长为单位）
K = 2; % 信号源个数

% 生成模拟接收信号数据（如有实际数据，请替换此部分）
N = 100; % 时间快拍数
angles = [-30, 20]; % 信号源的实际方向角（度）
lambda = 1; % 波长
theta_rad = angles * pi / 180; % 转换为弧度
A = exp(-1i * 2 * pi * d * (0:M-1).' * sin(theta_rad)); % 阵列流形矩阵
S = randn(K, N) + 1i * randn(K, N); % 信号源的信号
X = A * S + 0.1 * (randn(M, N) + 1i * randn(M, N)); % 接收信号加噪声

% 计算协方差矩阵
R = (X * X') / N;

% 特征值分解
[eigenvectors, eigenvalues] = eig(R);
[eigenvalues, idx] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:, idx);

% 噪声子空间
E_N = eigenvectors(:, K+1:end);

% 扫描角度范围
theta_scan = -90:0.1:90; % 扫描角度范围（度）
P_music = zeros(size(theta_scan)); % 初始化MUSIC谱

% 计算MUSIC谱
for i = 1:length(theta_scan)
    steering_vector = exp(-1i * 2 * pi * d * (0:M-1).' * sin(theta_scan(i) * pi / 180));
    P_music(i) = 1 / abs(steering_vector' * (E_N * E_N') * steering_vector);
end

% 找到谱峰
[~, peak_indices] = findpeaks(P_music);
[~, sorted_indices] = sort(P_music(peak_indices), 'descend');
doa_estimates = theta_scan(peak_indices(sorted_indices(1:K)));

% 绘制MUSIC谱
figure;
plot(theta_scan, 10*log10(P_music));
xlabel('Angle (degrees)');
ylabel('Spatial Spectrum (dB)');
title('MUSIC Spectrum');
grid on;

% 输出估计的方向角度
disp('Estimated DOA angles (degrees):');
disp(doa_estimates);
