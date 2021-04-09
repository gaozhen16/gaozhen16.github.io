function [W_rf, F_rf, H_eq] = Analog_SMSE2(H, Nt, Nr, Ns, K, Q, beta, rho)
% H_eq_k = W_rf_k' * H_k * F_rf;
%  DFT码本改为可过采样的，具体见Q_codebook
%  实现多个数据流传输;
%%
Mt = K * Ns; Mr = Ns; 

if length(Nr) == 2
    Nr_total = Nr(1) * Nr(2);
    Nt_total = Nt(1) * Nt(2);
else
    Nr_total = Nr;
    Nt_total = Nt;
end

DFT_set_r = Q_codebook(Nr, Q(1), rho); DFT_set_t = Q_codebook(Nt, Q(2), rho); % transmit and recieve codebook

len_r = size(DFT_set_r, 2); len_t = size(DFT_set_t, 2);
temp = zeros(len_r * K, len_t);
for i = 1 : K
    temp(((i - 1) * len_r + 1) : (i * len_r), :) ...
        = abs(DFT_set_r' * H(:, :, i) * DFT_set_t); % enery matrix
end

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
    
    dft_t_temp = abs(DFT_set_t(:, pos_t)' * DFT_set_t); % 基站去除相关性较强的码本
    dft_t_pos = dft_t_temp > beta;
    temp(:, dft_t_pos) = 0;
    
    dft_r_temp = abs(DFT_set_r(:, pos_rk)' * DFT_set_r); % 用户端去除相关性较强的码本
    dft_r_pos = find(dft_r_temp > beta);
    temp((k - 1) * len_r + dft_r_pos, :) = 0;
    
    if count_k(k) == 0
        error('天线数量较少(Nr, Nt)，导致码本相关性较强（a）');
    end
        
    W_rf(:, count_k(k), k) = DFT_set_r(:, pos_rk);
    F_rf(:, count_k(k), k) = DFT_set_t(:, pos_t); 
    if val_k >= Ns
        count_k(pos_k) = -1;
        temp(((pos_k - 1) * len_r + 1) : (pos_k * len_r), :) = 0;
    end
end

H_eq = zeros(Mr, Mt, K);
F_rf2 = reshape(F_rf, Nt_total, Mt);
for m = 1 : K
    H_eq(:, :, m) = W_rf(:, :, m)' * H(:, :, m) * F_rf2;
end