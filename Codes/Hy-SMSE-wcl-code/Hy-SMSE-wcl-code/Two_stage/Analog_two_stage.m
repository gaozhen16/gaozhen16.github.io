function [W_rf, F_rf, H_eq] = Analog_two_stage(H, Nt, Nr, Ns, K)
% 求出模拟预编码器和合并器
% H_eq_k = W_rf_k' * H_k * F_rf;
%%
Mt = K * Ns; Mr = Ns;
if length(Nr) == 2
    Nr_total = Nr(1) * Nr(2);
    Nt_total = Nt(1) * Nt(2);
else
    Nr_total = Nr;
    Nt_total = Nt;
end

% 生成码本
DFT_set_r = DFT_set(Nr); DFT_set_t = DFT_set(Nt); % transmit and recieve codebook

% 创建能量矩阵
len_r = size(DFT_set_r, 2); len_t = size(DFT_set_t, 2);
temp = zeros(len_r * K, len_t);
for i = 1 : K
    temp(((i - 1) * len_r + 1) : (i * len_r), :) ...
        = abs(DFT_set_r' * H(:, :, i) * DFT_set_t); % enery matrix
end

% 能量贪婪算法
W_rf = zeros(Nr_total, Mr, K); F_rf = zeros(Nt_total, Mr, K); 
count_k = zeros(1, K); %对应用户累计记录，当达到用户数据流上限，该用户清零
for i = 1 : Mt
    pos = find(temp == max(max(temp)));
    if length(pos) >= 2 % 出现多个最大值相同的情况
        pos = pos(1);
    end
    pos_r = mod(pos, K * len_r);% 确定用户端码本位置
    if pos_r == 0 % 当pos = K * len_r时
        pos_r = K * len_r;
    end
    k = ceil(pos_r / len_r); pos_rk = mod(pos_r, len_r);
    if pos_rk == 0 % 当pos_r = len_r时
        pos_rk = len_r;
    end
    pos_t = ceil(pos / (K * len_r));% 确定基站码本位置pos_t，用户k码本位置pos_rk

    temp(pos_r, :) = 0; % 对应位置清零，行、列、用户
    temp(:, pos_t) = 0;
    count_k(k) = count_k(k) + 1;  
    [val_k, pos_k] = max(count_k);


    W_rf(:, count_k(k), k) = DFT_set_r(:, pos_rk);
    F_rf(:, count_k(k), k) = DFT_set_t(:, pos_t); 
    if val_k >= Ns
        count_k(pos_k) = -1;
        temp(((pos_k - 1) * len_r + 1) : (pos_k * len_r), :) = 0;
    end
end

% 生成等效信道 H_eq
H_eq = zeros(Mr, Mt, K);
F_rf2 = reshape(F_rf, Nt_total, Mt);
for m = 1 : K
    H_eq(:, :, m) = W_rf(:, :, m)' * H(:, :, m) * F_rf2;
end