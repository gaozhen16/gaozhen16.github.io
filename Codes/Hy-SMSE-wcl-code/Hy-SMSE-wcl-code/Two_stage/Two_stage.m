function R = Two_stage(H, Ns, Nt, Nr, K, snr, msg, N_frame, noise, mod_obj, demod_obj)
% expression: y = W_RF ^H * (H * F_RF * F_BB * x + noise)
% 当 D 和 noise 没有输入时，R = 频谱效率；否则， R = 误码数量
%% set parameters  （only consider single stream）
if Ns > 1
    error('Ns > 1 for the Two-stage scheme')
end
Mr = 1;
Mt = K;
Ns = Mr;
P = Mt;
if length(Nt) == 2
    Nt_total = Nt(1) * Nt(2);
else
    Nt_total = Nt;
end
%% 
[W_RF, F_RF, H_eq] = Analog_two_stage(H, Nt, Nr, Ns, K); %
H_eq2 = zeros(K * Mr, Mt);
for m = 1 : K
    H_eq2(((m - 1) * Mr + 1) : (m * Mr), :) = H_eq(:, :, m);
end
F_BB = pinv(H_eq2);
% size(F_RF)
% size(F_BB)
F_RF2 = reshape(F_RF, Nt_total, Mt);
F_BB_a = sqrt(P) / norm(F_RF2 * F_BB,'fro');
F_BB = F_BB * F_BB_a;
%%
R = 0;
if nargin == 6
    for m = 1 : K
        signal = (P / K / Ns) * W_RF(:, :, m)' * H(:, :, m) * F_RF2 * F_BB(:, m) * F_BB(:, m)' * F_RF2' * H(:, :, m)' * W_RF(:, :, m);
        noise =(P / snr) * W_RF(:, :, m)' * W_RF(:, :, m);
        interference = 0;
        for m_i = 1 : K
            if m_i ~= m
                interference = interference + W_RF(:, :, m)' * H(:, :, m) * F_RF2 * F_BB(:, m_i) * F_BB(:, m_i)' * F_RF2' * H(:, :, m)' * W_RF(:, :, m);
            end
        end
        noise = noise + interference;
        R = R + log2(abs(det(eye(Ns) + signal * pinv(noise))));
    end
else
    symbol = modulate(mod_obj, msg);
    Scale = modnorm(symbol,'avpow',1);
    Symbol_nomalized = reshape(Scale * symbol, K * Ns, N_frame); % 归一化
    y = [];
    for m = 1 : K
        y = [y; W_RF(:, :, m)' * (H(:, :, m) * F_RF2 * F_BB * Symbol_nomalized + noise)];
    end
    Symbol_hat = reshape(y / Scale / F_BB_a, K * Ns * N_frame, 1);  %当存在幅度上的调试时，需要考虑幅度
    msg_hat = demodulate(demod_obj, Symbol_hat);
    R = sum(msg_hat ~= msg);
end