% code: CS-based Active User Detection and Channel Estimation for Massive Access (DSAMP)
% spatial-frequency structured sparsity
% written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: 2019.03.25
% sim time:  N_sim = 3000, time = 4h, 
% predicted: N_sim = 3000, time = 4h, 1903  -- 1903

%% Prepare the workspace && Set parameters 
clc; clear; close all
% parameters of system && channel model
K = 500;           % number of potential UDs in the cell
Ka = 50;           % number of active UDs
N_BS = 64;         % number of BS antennas
P = 1;             % number of pilot subcarriers
SNR = 30;          % SNR in dB
T_max = K;         % maximum time slot overhead
N = 2048;          % number of subcarriers
Bw = 10e6;         % system bandwidth
fs = Bw;           % frequency of sampling
L_cp = 64;         % length of cyclic prefix
sigma2_alpha = 1;  % average power of path, generally, sigma2_alpha = 1;
Lp_min = 8;        % number of paths
Lp_max = 14;
% parameters of AUD
p_cg = 0.9;

%% Simulations
tic
N_sim = 3000;
T_set = 30:2:80;                   % set time slot overhead
Pe_cg = zeros(numel(T_set),1);     % EDP of CG-AD
MSE_dsamp = zeros(numel(T_set),1);  
startmatlabpool(26)
for sim_ii = 1:N_sim
   %% Generate simulation samples
    % generate channel model
%     H = zeros(K,N_BS,P);
%     index = randperm(K);
%     H(index(1:Ka),:,:) = sqrt(1/2).*(randn(Ka,N_BS,P)+1i.*randn(Ka,N_BS,P));
    act_flag = zeros(K,1);       
    index = randperm(K);
    act_flag(index(1:Ka)) = 1; 
    H = zeros(K,N_BS,P);
    for k = 1:K 
        Lp = randi([Lp_min,Lp_max],1,1);
        [H_f, H_ang_f] = generate_FSF_channel_grid_match(N_BS, N, P, Bw, fs, L_cp, sigma2_alpha, Lp);
        H(k,:,:) = act_flag(k).*H_f;
    end
    % pilot matrix 
    gamma = 1;
    S = sqrt(gamma/2).*(randn(T_max,K,P)+1i.*randn(T_max,K,P));
    % noisy received signals
    Y = zeros(T_max,N_BS,P);
    for p = 1:P
        Y(:,:,p) = awgn(S(:,:,p)*H(:,:,p), SNR, 'measured');
    end
     
   %% Active User Detection and Channel Estimation
    for T_ii = 1:numel(T_set)
        T = T_set(T_ii);
        
        % DMMV-AMP
        TH = 0.005;
        [H_dsamp, supp, iter] = SAMP_distributed(Y(1:T,:,:), S(1:T,:,:), TH);
        % CG-AD
        act_cg = zeros(K,1);
        epsilon_cg = 0.01*max(max(max(abs(H_dsamp))));
        act_cg(sum(sum(abs(H_dsamp)>epsilon_cg,3),2)./N_BS./P >= p_cg) = 1;        
        
        % Data recording
        Pe_cg(T_ii) = Pe_cg(T_ii) + sum(abs(act_cg-act_flag))/K;
        MSE_dsamp(T_ii) = MSE_dsamp(T_ii) + norm(H_dsamp(1:end)-H(1:end))^2/(K*N_BS*P);
        
        fprintf('sim = %d, T = %d, Pe_cg = %4.5f, MSE_sp = %7.9f\n', sim_ii, T, ...
                sum(abs(act_cg-act_flag))/K, norm(H_dsamp(1:end)-H(1:end))^2/(K*N_BS*P));
        
    end 
end
closematlabpool
toc

save Data_DSAMP_M64_P1.mat K Ka N_BS P SNR Lp_min ...
     Lp_max T_set Pe_cg MSE_dsamp N_sim

Pe_cg = Pe_cg./N_sim;
MSE_dsamp = MSE_dsamp./N_sim;

% figure Pe
figure
semilogy(T_set, Pe_cg, 'b-x', 'linewidth', 1.5);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead');ylabel('P_e');
legend('CG-AD');

% figure MSE
figure
semilogy(T_set, MSE_dsamp, 'r--o', 'linewidth', 1.5);
hold on
set(gca, 'GridLineStyle', '-.');
grid on;
xlabel('Time Overhead');ylabel('MSE');
legend('DMMV-OMP');





















