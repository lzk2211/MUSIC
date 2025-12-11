N = 4;    % 采样天线数量，空间轴维度
M = 2;    % 空间上有多少个信号源
T = 1024; % 采样时间点，时间轴维度
d = 0.06; % 天线间距

signal_F = 2.4e9;   % 发射源信号的频率
lambda = 3e8 / signal_F; % 0.125m

angles = [-10, 50]; % 信号源的实际方向角（度）
theta_rad = angles * pi / 180; % 转换为弧度

% 信号向量exp(1i*(0:N-1).'*Ω)
% Ω = 2*pi*d/lambda*sin(theta_rad)代表空间上的频率，即信号到达天线的时刻不同导致相位的变化频率
A = exp(-1i * 2 * pi * d/lambda * (0:N-1).' * sin(theta_rad)); % 信号向量，天线数量N * 频率源数量M

f = [20, 10]; % 正弦波频率，多个目标的话，不能是同频率否则就混叠一起了
t = (0:T-1).'/T; % 时间戳，T其实是fs
S = exp(1i * 2*pi * t * f); % 信号波形，对于MUSIC算法可以是任意信号

% X = A * S.' + randn(N, T); % 4*2 * 2*100 得到N*T的采样数据 + 噪声
X = awgn(A * S.', 0); % 4*2 * 2*100 得到N*T的采样数据 + 噪声

% 计算协方差矩阵
R = (X * X') / T;          % X信号进行自相关，多个采样时间T下的信号做完自相关后求平均值

% 特征值分解
[eigenvectors, eigenvalues] = eig(R); % 特征值中，有M个值会明显比N-M个值大很多
[eigenvalues, idx] = sort(diag(eigenvalues), 'descend'); % 将大的值排序到前面
eigenvectors = eigenvectors(:, idx); % 得到对应排序的特征向量矩阵，前面M是信号子空间，后面N-M是噪声子空间

% 噪声子空间
E_N = eigenvectors(:, M+1:end); % 求噪声子空间为后续MUSIC做准备

% 扫描角度范围
theta_scan = -90:0.1:90; % 扫描角度范围（度）
P_music = zeros(size(theta_scan)); % 初始化MUSIC谱

% 计算MUSIC谱
for i = 1:length(theta_scan) % 对附近的角度的可能空间频率信号的旋转因子进行遍历匹配
    steering_vector = exp(-1i * 2 * pi * d/lambda * (0:N-1).' * sin(theta_scan(i) * pi / 180));
    % abs(steering_vector' * (E_N * E_N') * steering_vector)
    % 扫描的信号旋转因子向量，在接收到的信号的噪声子空间中的投影向量和该信号旋转因子向量的自相关函数
    % 如果没有噪声子空间上的投影，说明扫描到的信号是信号子空间中的信号，那么该空间扫描频率是信号源的空间频率，该式子为非常小的值
    P_music(i) = 1 / abs(steering_vector' * (E_N * E_N') * steering_vector);
end

% 找到谱峰
[~, peak_indices] = findpeaks(P_music); % 寻找扫描结果中的最大值
[~, sorted_indices] = sort(P_music(peak_indices), 'descend');
doa_estimates = theta_scan(peak_indices(sorted_indices(1:M))); % 在所有可能的结果中选择前M个结果认定为最后结果

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
