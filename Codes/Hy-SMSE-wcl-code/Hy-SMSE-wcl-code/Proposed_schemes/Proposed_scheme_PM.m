function R = Proposed_scheme_PM(H, Nt, Nr, Ns, K, Q, snr, msg, N_frame, noise, mod_obj, demod_obj)
% y_k = W_bb_k' * W_rf_k' * (H * F_rf * F_bb * x + n_k)
% Proposed Digital Part Using Analog Part in [8]
%%
W_rf = RF_combiner(H, Nr, Ns, K, Q(1));
[F_rf, H_eq] = RF_precoder(H, W_rf, K, Ns, Nt, Q(2));

V_BD = channl_SVD(H_eq, K, Ns);
[F_bb, W_bb, W_a] = Digital_SMSE(Nt, K, Ns, H_eq, V_BD, snr, W_rf, F_rf);
S = zeros(1, K * Ns) + W_a;
%%
if nargin <= 7
    R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr);
else
    R = Ber(Ns, K, W_bb, W_rf, F_bb, F_rf, H, S, msg, N_frame, noise, mod_obj, demod_obj);
end