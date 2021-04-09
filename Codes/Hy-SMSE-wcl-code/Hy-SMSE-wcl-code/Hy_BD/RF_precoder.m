function [F_rf, H_eq] = RF_precoder(H, W_rf, K, Ns, Nt, Bt)
% create the RF precoder at BS
%% create the effective channel H_int
Mt = K * Ns; Mr = Ns; 
if length(Nt) == 2
    Nt_total = Nt(1) * Nt(2);
else
    Nt_total = Nt;
end
H_int = zeros(Nt_total, Mr, K); % -> 由于reshape函数按列抽取，需要转置
for m = 1 : K
    M_k = W_rf(:, :, m);
    H_k = H(:, :, m);
    H_int(:, :, m) = (M_k' * H_k).';
end
H_int2 = reshape(H_int, Nt_total, K * Mr).'; % 将三维数组转化为二维

%% create and quantify F
if nargin == 5
    F_rf = exp(1i * angle(H_int2')) / sqrt(Nt_total);
else
    B = 2 ^ Bt;
    F_rf = round((angle(H_int2) + pi) * B / 2 / pi) * 2 * pi / B;
    F_rf = (exp(1i * F_rf) / sqrt(Nt_total))';
end

%% create the equivalent channel H_eq
H_eq = zeros(Mr, Mt, K);
for i = 1 : K
    H_eq(:, :, i) = H_int(:, :, i).' * F_rf;
end