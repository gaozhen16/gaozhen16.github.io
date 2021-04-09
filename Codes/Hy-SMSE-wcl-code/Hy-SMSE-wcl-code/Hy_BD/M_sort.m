function [M_k, norm1] = M_sort(H_k, D, D_len, Mr)
% 排序求出每个用户的W矩阵
% D 码本
% D_len 码本大小
% norm1 存储每个d * H的范数和序号(按phase排序的)
%%
dH_norm = zeros(2,D_len); % 存储每个d * H的范数和序号
for i = 1 : D_len
    A = D(:, i)' * H_k;   
    dH_norm(1, i) = norm(A,'fro');
end
dH_norm(2, :) = (1 : D_len);
norm1 = dH_norm(1, :);

% 排序(冒泡排序, 降序)
for i = 1 : D_len
    for j = 1 : (D_len - i)
        if dH_norm(1,j) < dH_norm(1,j + 1)
            a = dH_norm(:,j);
            dH_norm(:,j) = dH_norm(:,j + 1);
            dH_norm(:,j + 1) = a;
        end
    end
end
M_k = zeros(length(D(:, 1)), Mr);
for i = 1 : Mr
    M_k(:, i) = D(:,dH_norm(2, i));
end
