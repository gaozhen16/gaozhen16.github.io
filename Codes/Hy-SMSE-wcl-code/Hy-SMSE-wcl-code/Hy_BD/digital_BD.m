function [M, B, S, B_a] = digital_BD(H_eq, F, K, Ns)
% 利用等效信道两次SVD分解，求出B、M
%%
Mt = K * Ns; Mr = Ns;Pt = Mt;
M = zeros(Mr, Ns, K);
B = zeros(Mt, Ns, K);
S = zeros(1, K * Ns);
H_eq2 = zeros(K * Mr, Mt);
for m = 1 : K
    H_eq2(((m - 1) * Mr + 1) : (m * Mr), :) = H_eq(:, :, m);
end
for m = 1 : K
    % define H_
    if m == 1
        H_ = H_eq2(Mr * m + 1 : end, :);
    else
        H_1 = H_eq2(1 : Mr * (m - 1), :);
        H_2 = H_eq2(Mr * m + 1 : end, :);
        H_ = [H_1 ; H_2];
    end
 
    % 两次SVD分解
    [~,~,V1] = svd(H_);
    V1_ = V1(:,K * Mr - Mr + 1 : end);
    H_k = H_eq2(m * Mr - Mr + 1 : m * Mr, :);  %到k用户的等效信道
    [U1, S1, V1] = svd(H_k * V1_);
    % obtain B and M
    S(1, ((m - 1) * Ns + 1) : (m * Ns)) = diag(S1);
    V1 = V1(:,1 : Ns);
    U1 = U1(:,1 : Ns);
    M(:, :, m) = U1;
    B(:, :, m) = V1_ * V1;  
end
B2 = reshape(B, Mt, K * Ns);
if size(F) == [1 1]  % 临时处理
    B_a = sqrt(Pt) / norm(B2, 'fro' );
else
    B_a = sqrt(Pt) / norm(F * B2, 'fro' );
end
B2 = B2 * B_a;  % 功率归一化处理
B = reshape(B2, Mt, Ns, K);
S = S * B_a;
