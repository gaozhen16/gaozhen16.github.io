function [F_bb, W_bb, W_a] = Digital_SMSE(Nt, K, Ns, H_eq, V_BD, snr, W_rf, F_rf)
% 由闭合表达式求出F_bb, W_bb
%%
Mt = K * Ns; Mr = Ns; Pt = Mt;
if length(Nt) == 2
    Nt_total = Nt(1) * Nt(2);
else
    Nt_total = Nt;
end

F_rf2 = reshape(F_rf, Nt_total, Mt);
H_eq2 = zeros(K * Mr, Mt);
for m = 1 : K
    H_eq2(((m - 1) * Mr + 1) : (m * Mr), :) = H_eq(:, :, m);
end
[~, S, ~] = svd(H_eq2);
% if min(diag(S)) < 1
%     H_eq2
%     S
% % elseif min(diag(S)) < 2
% %     S
% end

F_bb2 = pinv(H_eq2' * V_BD' * V_BD * H_eq2) * H_eq2' * V_BD';
W_a = sqrt(Pt) / norm(F_rf2 * F_bb2, 'fro');
% xi = trace(V_BD' * H_eq2 * F_bb2 * F_bb2' * H_eq2' * V_BD) - trace(V_BD' * H_eq2 * F_bb2) - trace(F_bb2' * H_eq2' * V_BD);
% abs(xi)

% xi = trace(V_BD' * H_eq2 * F_bb2 * F_bb2' * H_eq2' * V_BD) - trace(V_BD' * H_eq2 * F_bb2) - trace(F_bb2' * H_eq2' * V_BD);
% abs(xi)
% mm
F_bb = reshape(F_bb2, Mt, Ns, K); % ->转化为3维数组
% if isWa == 1 % 当is_W_a = 1时
%     W_a = 1; 
% end
    
W_bb = zeros(Ns, Mr, K);
for m = 1 : K
    H_eq_k = H_eq(:, :, m);
    F_bb_k = F_bb(:, :, m);
    M_rf_k = W_rf(:, :, m);
    W_bb(:, :, m) = (F_bb_k' * H_eq_k' * pinv(H_eq_k * F_bb2 * F_bb2' * ...
        H_eq_k' + Pt / snr / (W_a^2) * M_rf_k' * M_rf_k))';
end
F_bb = F_bb * W_a;   % 功率归一化

