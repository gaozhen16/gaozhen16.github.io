function R = Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, beta, rho, msg, N_frame, noise, mod_obj, demod_obj)
% y_k = W_bb_k' * W_rf_k' * (H * F_rf * F_bb * x + n_k)
% 基于（DFT+DFT）能量最大方案
% DFT码本改为可过采样的，具体见Q_codebook
% 测试多流传输
%%
[W_rf, F_rf, H_eq] = Analog_SMSE2(H, Nt, Nr, Ns, K, Q, beta, rho);

% V_BD = channl_SVD(H_eq, K, Ns);
V_BD = eye(K * Ns);
[F_bb, W_bb, W_a] = Digital_SMSE(Nt, K, Ns, H_eq, V_BD, snr, W_rf, F_rf);
S = zeros(1, K * Ns) + W_a;
% W_bb(:, :, 1)' * H_eq(:, :, 1) * F_bb(:, :, 1)
if nargin <= 10
    R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
else
    R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, S, msg, N_frame, noise, mod_obj, demod_obj);
end