function W_rf = RF_combiner(H, Nr, Ns, K, Br)
% create the RF combiners at the users
%% create user codebook
if nargin == 4
    D = DFT_set(Nr);
else
    D = DFT_set(Nr, Br);
end

%% create M
if length(Nr) == 2
    Nr_total = Nr(1) * Nr(2);
else
    Nr_total = Nr;
end
W_rf = zeros(Nr_total, Ns, K);  % define M's size
% norm1 = zeros(1, length((D(1, :))), K);
for m = 1 : K
    H_k = H(:, :, m);
    [W_rf(:, :, m), ~] = M_sort(H_k, D, length(D(1, :)), Ns);
end