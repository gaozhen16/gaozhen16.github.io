function R = Hy_BD(H, Nt, Nr, Ns, K, Q, snr, msg, N_frame, noise, mod_obj, demod_obj)
% BD 混合预编码
% expression: y = {W_bb}^H * {W_rf}^H * (H * F_rf * F_bb * x + noise)
% 当 D 和 noise 没有输入时，R = 频谱效率；否则， R = 误码数量
%%
W_rf = RF_combiner(H, Nr, Ns, K, Q(1));
[F_rf, H_eq] = RF_precoder(H, W_rf, K, Ns, Nt, Q(2));

% digital precoding
[W_bb, F_bb, S, ~] = digital_BD(H_eq, F_rf, K, Ns);

%%
if nargin <= 7
    R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
else
    R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, S, msg, N_frame, noise, mod_obj, demod_obj);
end
