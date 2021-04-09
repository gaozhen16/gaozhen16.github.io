function R = spectral_efficiency(Ns, K, W_bb, W_rf, F_bb, H_eq, snr)
% 计算总的频谱效率
%%
Pt = K * Ns; R = 0;
for m = 1 : K
    W_bb_k = W_bb(:, :, m);
    H_eq_k = H_eq(:, :, m);
    F_bb_K = F_bb(:, :, m);
    W_rf_k = W_rf(:, :, m);

    VH = W_bb_k' * H_eq_k;
    signal = (Pt / K / Ns) * VH * F_bb_K * F_bb_K' * VH';
    noise = (Pt / snr) * W_bb_k' * W_rf_k' * W_rf_k * W_bb_k;
    interference = 0;
    for m_i = 1 : K
        if m_i ~= m
            F_bb_ki = F_bb(:, :, m_i);
            interference = interference + VH * F_bb_ki * F_bb_ki' * VH';
        end
    end
    interference = (Pt / K / Ns) * interference;
    noise = noise + interference;
%     fprintf('\nHy-SMSE %f %f\n', abs(det(signal)), abs(det(noise)));
    R = R + log2(abs(det(eye(Ns) + signal * pinv(noise))));       
end