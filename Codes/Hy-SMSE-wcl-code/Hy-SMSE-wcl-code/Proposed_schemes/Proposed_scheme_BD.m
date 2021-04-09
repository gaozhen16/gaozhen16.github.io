function R = Proposed_scheme_BD(H, Nt, Nr, Ns, K, snr, Q, beta, rho, msg, N_frame, noise, mod_obj, demod_obj)
% y_k = W_bb_k' * W_rf_k' * (H * F_rf * F_bb * x + n_k),数字BD
% DFT码本改为可过采样的，具体见Q_codebook
% 测试多流传输
%%
Mt = K * Ns;
if length(Nt) == 2
    Nt_total = Nt(1) * Nt(2);
else
    Nt_total = Nt;
end

[W_rf, F_rf, H_eq] = Analog_SMSE2(H, Nt, Nr, Ns, K, Q, beta, rho);
F_rf2 = reshape(F_rf, Nt_total, Mt);
[W_bb, F_bb, S, ~] = digital_BD(H_eq, F_rf2, K, Ns);

if nargin <= 9
    R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
else
    R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, S, msg, N_frame, noise, mod_obj, demod_obj);
end