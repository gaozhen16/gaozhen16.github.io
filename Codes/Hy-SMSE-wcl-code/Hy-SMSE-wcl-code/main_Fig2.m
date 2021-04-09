clc; clear; add_path();

Nt = [8 8]; Nr = [4 4]; K = 8; Ns = 1;
Q = [log2(Nr(1)) log2(Nt(1))];
SNR = (-10 : 4 : 10); 
channel_N = 500; % 实验次数 

R = zeros(13, length(SNR));
for i = 1 : length(SNR)
    snr =10.^(SNR(i) / 10);
    for channel_i = 1 : channel_N   
        tic
        
        [H, H_At, H_Ar] = channel(Nt, Nr, K);
        R(1, i) = R(1, i) + Two_stage(H, Ns, Nt, Nr, K, snr); 
        R(2, i) = R(2, i) + Hy_BD(H, Nt, Nr, Ns, K, Q, snr); 
        R(3, i) = R(3, i) + Proposed_scheme_PM(H, Nt, Nr, Ns, K, Q, snr); 
        R(4, i) = R(4, i) + Proposed_scheme_BD(H, Nt, Nr, Ns, K, snr, Q, 0.15, 8); 
        R(5, i) = R(5, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 1, 1); 
        R(6, i) = R(6, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 1, 2); 
        R(7, i) = R(7, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 1, 8); 
        R(8, i) = R(8, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, 2); 
        R(9, i) = R(9, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, 4); 
        R(10, i) = R(10, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, 8); 
        R(11, i) = R(11, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, -8); % \rho < 0，表示无量化
        R(12, i) = R(12, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, -16);
        toc
    end
    %----------------------------------完成度
    disp([num2str(i * 100 / length(SNR)),'%']);
    %----------------------------------------------
end
R = R / channel_N;

l_width = 1.5;

plot(SNR, R(1, :), 's-m', 'LineWidth', l_width);
grid on, hold on
plot(SNR, R(2, :), '+-g', 'LineWidth', l_width);
plot(SNR, R(3, :), '*-r', 'LineWidth', l_width);
plot(SNR, R(4, :), 'd-k', 'LineWidth', l_width);
plot(SNR, R(5, :), '+-.b', 'LineWidth', l_width);
plot(SNR, R(6, :), 'v-b', 'LineWidth', l_width);
plot(SNR, R(7, :), '*-b', 'LineWidth', l_width);
plot(SNR, R(8, :), '+-b', 'LineWidth', l_width);
plot(SNR, R(9, :), 'o-b', 'LineWidth', l_width);
plot(SNR, R(10, :), 's-b', 'LineWidth', l_width);
plot(SNR, R(11, :), '-.c', 'LineWidth', l_width);
plot(SNR, R(12, :), '-.b', 'LineWidth', l_width);


legend('Two-Stage Scheme [6]', 'Hybrid BD Scheme [8]', ...
    'Proposed Digital Part Using Analog Part in [8]', ...
    'Proposed Analog Part (\rho = 8,\beta = 0.15) Using Digital Part in [8]',...
    'Proposed Scheme, \rho = 1, 0 < \beta \leq 1', ...
    'Proposed Scheme, \rho = 2, \beta = 1', ...
    'Proposed Scheme, \rho = 8, \beta = 1', ...
    'Proposed Scheme, \rho = 2, \beta = 0.15',...
    'Proposed Scheme, \rho = 4, \beta = 0.15', ...
    'Proposed Scheme, \rho = 8, \beta = 0.15', ...
    'Proposed Scheme with B_t = B_r = \infty, \rho = 8,\beta = 0.15', ...
    'Proposed Scheme with B_t = B_r = \infty, \rho = 16,\beta = 0.15')
xlabel('SNR(dB)'),ylabel('SSE');
set(gcf, 'position', [650 250 600 460])
set(gca, 'xtick', SNR);
% save('Fig2_180406_3')

