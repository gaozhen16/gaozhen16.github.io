function [R, error] = opt_scheme(H, H_At, H_Ar, Nt, Nr, Ns, K, snr, choose, beta, msg, N_frame, noise, mod_obj, demod_obj)
% infinite resolution
% clc; clear; add_path();
% Nt = 16; Nr = 16; K = 2; Ns = 2;
error = 0;  % 防止错误
if length(Nt) == 2
    Nt_total = Nt(1) * Nt(2);
    Nr_total = Nr(1) * Nr(2);
else
    Nt_total = Nt;
    Nr_total = Nr;
end

% [H, H_At, H_Ar] = channel(Nt, Nr, K);

% Analog part
Mt = K * Ns; Mr = Ns;
W_rf = zeros(Nr_total, Mr, K);
F_rf = zeros(Nt_total, Mr, K);

len_r = size(H_Ar(:, :, 1), 2); 
% len_t = size(H_At(:, :, 1), 2);
temp = zeros(1, K * len_r);
H_Ar_all = zeros(Nr_total, K * len_r);
H_At_all = zeros(Nt_total, K * len_r);
for m = 1 : K
    H_Ar_all(:, ((m - 1) * len_r + 1) : (m * len_r)) = H_Ar(:, :, m);
    H_At_all(:, ((m - 1) * len_r + 1) : (m * len_r)) = H_At(:, :, m);
    temp(((m - 1) * len_r + 1) : (m * len_r)) = diag(abs(H_Ar(:, :, m)' * H(:, :, m) * H_At(:, :, m)));
end
% beta = 0.5; % 设置参数，降低相关性
count_k = zeros(1, K); %对应用户累计记录，当达到用户数据流上限，该用户清零
for i = 1 : (K * Ns)
    pos = find(temp == max(temp));
    if length(pos) > 1  % 不完善的地方，当temp为0向量时，可能存在误判
        pos = pos(1);
    end
    k = ceil(pos / len_r); % 向上取整
    count_k(k) = count_k(k) + 1;
    
    if count_k(k) == 0
        error = 1;
        break;
        error('--> /beta');
    end

    W_rf(:, count_k(k), k) = H_Ar_all(:, pos);
    F_rf(:, count_k(k), k) = H_At_all(:, pos);
    
    H_At_temp = abs(H_At_all(:, pos)' * H_At_all); % 去除相关性较强的，基站端
    H_At_pos = H_At_temp > beta;
    temp(H_At_pos) = 0;
    
    if count_k(k) >= Ns
        count_k(k) = -1;
        temp(((k - 1) * len_r + 1) : (k * len_r)) = 0;
    else
        H_Ar_temp = abs(H_Ar_all(:, pos)' * H_Ar_all(:, ((k - 1) * len_r + 1) : (k * len_r)));
        H_Ar_pos = find(H_Ar_temp > beta);
        temp(H_Ar_pos + (k - 1) * len_r) = 0;
    end
    
    
end


    % Digital part
    H_eq = zeros(Mr, Mt, K);
    % size(F_rf)
    F_rf2 = reshape(F_rf, Nt_total, Mt);
    for m = 1 : K
        H_eq(:, :, m) = W_rf(:, :, m)' * H(:, :, m) * F_rf2;
    end
    H_eq2 = zeros(K * Mr, Mt);
    for m = 1 : K
        H_eq2(((m - 1) * Mr + 1) : (m * Mr), :) = H_eq(:, :, m);
    end
    [~, S, ~] = svd(H_eq2);
    if (error ~= 1) && (min(diag(S)) < 1)
        error = 1;
        disp('min < 1 ----->');
    end
if error == 0
    if choose == 1 %不同的数字方案
        V_BD = channl_SVD(H_eq, K, Ns);
        [F_bb, W_bb, W_a] = Digital_SMSE(Nt, K, Ns, H_eq, V_BD, snr, W_rf, F_rf);
        S = zeros(1, K * Ns) + W_a;
    elseif choose == 2
        [W_bb, F_bb, S, ~] = digital_BD(H_eq, F_rf2, K, Ns);
    end

    if nargin <= 10
        R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
    else
        R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, S, msg, N_frame, noise, mod_obj, demod_obj);
    end
else
    R = 0;
end
% if nargin <= 9
%     R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
% else
%     R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, msg, N_frame, noise, mod_obj, demod_obj);
% end