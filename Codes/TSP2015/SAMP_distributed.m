function [H_hat, support, support_num] = SAMP_distributed(Y, Phi, TH)
% SAMP: Sparsity Adaptive Matching Pursuit algoritm for compressed sensing.
% For theoretical analysis, please refer to the paper : Thong. T. Do, Lu Gan and Trac D. Tran ,、
% "Sparsity Adaptive Matching Purusit for practical compressed sensing" available at http://dsp.ece.rice.edu/cs

% Written by Thong Do(thongdo@jhu.edu)
% Updated on July, 26th 2008

% stop 条件改变为某一Block联合MMV分量的norm小于一定门限

% parameter usage:
%   Y: received signals 
%   Phi: measurement matrix
%   H_hat: reconstructed sparse signal
%   iter_num: number of iterations

[G, M, P] = size(Y); 
K = size(Phi,2);
%% Initialization
T_cur_stage = 1; % sparsity of current stage
i = 1;           % iteration index
j = 1;           % stage index
supp_set = [];
b_res = Y;
energy_b_last = inf;
iter_num = 0;
while  1 
    b_res_pre = b_res;
    % candidate list
    a_inn_pro = zeros(K,M,P);
    for pp = 1:P 
        a_inn_pro(:,:,pp) = Phi(:,:,pp)'*b_res(:,:,pp);
    end
    [~, sort_index] = sort(sum(sum(abs(a_inn_pro).^2,3),2), 'descend');
    candidate_set = union(supp_set, sort_index(1:T_cur_stage));
    % LS
    t_ls = zeros(K,M,P);
    for pp = 1:P
        t_ls(candidate_set,:,pp) = Phi(:,candidate_set,pp) \ Y(:,:,pp);
    end
    % purning support
    [~, sort_index2] = sort(sum(sum(abs(t_ls).^2,3),2), 'descend');
    Omiga = sort_index2(1:T_cur_stage);
    % LS new and residual 
    c_ls = zeros(K,M,P);
    for pp = 1:P
        c_ls(Omiga,:,pp) = Phi(:,Omiga,pp) \ Y(:,:,pp);
        b_res(:,:,pp) = Y(:,:,pp) - Phi(:,:,pp)*c_ls(:,:,pp);
    end
    [~, index_min] = min(sum(sum(abs(c_ls(Omiga,:,:)).^2,3),2));
    L_min = Omiga(index_min);
    if sum(sum(abs(c_ls(L_min,:,:)).^2,3),2)/(P*M) < TH
       break;
    elseif energy_b_last < sum(sum(sum(abs(b_res).^2,3),2))
       break;
    elseif sum(sum(sum(abs(b_res_pre).^2,3),2)) <= sum(sum(sum(abs(b_res).^2,3),2))
        j = j + 1;
        T_cur_stage = j;
        b_last = b_res;
        energy_b_last = sum(sum(sum(abs(b_last).^2,3),2));
    else
        supp_set = Omiga;
        i = i + 1;
    end
    iter_num = iter_num +1;  
end 
 
% reconstruction
H_hat = zeros(K,M,P);
for pp = 1:P
    H_hat(supp_set,:,pp) = Phi(:,supp_set,pp)\Y(:,:,pp);
end
support = sort(supp_set,'ascend');
support_num = numel(support);
end