clc; clear; warning off; add_path(); 
Nt = [8 8]; Nr = [4 4]; K = 2; Ns = 2; 
Pt = K * Ns; Q = [log2(Nr(1)) log2(Nt(1))];
SNR = (-12 : 4 : 4);
if length(Nr) == 2
    Nr_total = Nr(1) * Nr(2);
else
    Nr_total = Nr;
end

b_s = 4; % number of bits per symbol (1: 2QAM; 2: 4QAM; ...)
N_frame = 60;                   %帧数
N_packet = 10000;                 %分组数
error_N = 600;       % ber 累计错误上限
Max_ii = 10;
M = 2 ^ b_s;              % modulation order
N_pbits = N_frame * K * Ns * b_s;
N_tbits = N_pbits * N_packet;

mod_obj = modem.qammod('M',M,'SymbolOrder','Gray','InputType','bit');%调制方式设置
demod_obj = modem.qamdemod(mod_obj);

R = zeros(10, length(SNR));E = 0; opt = 0.35;
for i = 1 : length(SNR)
    snr =10.^(SNR(i) / 10);
    rng(0,'twister');
%     tic
    for N_packet_i = 1 : N_packet
        tic
        msg = randi([0,1],[N_pbits, 1]);  % message
        noise = (randn(Nr_total, N_frame) + 1i * randn(Nr_total, N_frame)) / sqrt(2);
        noise = sqrt(Pt / snr) * noise;
        for ii = 1 : Max_ii % 曲线更光滑
            [H, H_At, H_Ar] = channel(Nt, Nr, K);
            [temp_R7, error] = opt_scheme(H, H_At, H_Ar, Nt, Nr, Ns, K, snr, 1, opt, msg, N_frame, noise, mod_obj, demod_obj);
            if error == 1
                E = E + 1; %fprintf('error = %d\n', E);
            else
                break;
            end
        end
        R(1, i) = R(1, i) + Hy_BD(H, Nt, Nr, Ns, K, Q, snr, msg, N_frame, noise, mod_obj, demod_obj);
        R(2, i) = R(2, i) + Proposed_scheme_PM(H, Nt, Nr, Ns, K, Q, snr, msg, N_frame, noise, mod_obj, demod_obj);
        R(3, i) = R(3, i) + Proposed_scheme_BD(H, Nt, Nr, Ns, K, snr, Q, 0.15, 8, msg, N_frame, noise, mod_obj, demod_obj);
        R(4, i) = R(4, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, 8, msg, N_frame, noise, mod_obj, demod_obj);
%         R(5, i) = R(5, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.2, 8, msg, N_frame, noise, mod_obj, demod_obj);
        R(6, i) = R(6, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.3, 8, msg, N_frame, noise, mod_obj, demod_obj);
        R(7, i) = R(7, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.4, 8, msg, N_frame, noise, mod_obj, demod_obj);
%         R(7, i) = R(7, i) + temp_R7; 
%         R(8, i) = R(8, i) + opt_scheme(H, H_At, H_Ar, Nt, Nr, Ns, K, snr, 2, opt, msg, N_frame, noise, mod_obj, demod_obj); 
        R(9, i) = R(9, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, -8, msg, N_frame, noise, mod_obj, demod_obj); % \rho < 0，表示无量化
        R(10, i) = R(10, i) + Proposed_scheme(H, Nt, Nr, Ns, K, snr, Q, 0.15, -16, msg, N_frame, noise, mod_obj, demod_obj);
%         toc
        if R(10, i) >= error_N
            R(:, i) = R(:, i) / N_pbits / N_packet_i;
            break;
        elseif N_packet_i >= N_packet
             R(:, i) = R(:, i) / N_pbits / N_packet_i;
        end
    end   
%     save('Fig_4')
    %----------------------------------完成度
    disp([num2str(i * 100 / length(SNR)),'%, iter = ' num2str(N_packet_i)]);
    %----------------------------------------------
end

l_width = 1.5;
semilogy(SNR, R(1, :),'+-g','LineWidth',l_width)
grid on
hold on
semilogy(SNR, R(2, :),'s-r','LineWidth',l_width)
semilogy(SNR, R(3, :),'d-k','LineWidth',l_width)
semilogy(SNR, R(4, :),'o-b','LineWidth',l_width)
% semilogy(SNR, R(5, :),'+-b','LineWidth',l_width)
semilogy(SNR, R(6, :),'*-b','LineWidth',l_width)
semilogy(SNR, R(7, :),'s-b','LineWidth',l_width)%     
% semilogy(SNR, R(8, :),'-.r','LineWidth',l_width)
semilogy(SNR, R(9, :),'-.c','LineWidth',l_width)
semilogy(SNR, R(10, :),'-.b','LineWidth',l_width)
xlabel('SNR(dB)'),ylabel('BER');
legend('Hybrid BD Scheme [8]', ...
    'Proposed Digital Part Using Analog Part in [8]',...
    'Proposed Analog Part (\rho = 8, \beta = 0.2) Using Digital Part in [8]',...
    'Proposed Scheme, \rho = 8, \beta = 0.15', ...
    'Proposed Scheme, \rho = 8, \beta = 0.3', ...
    'Proposed Scheme, \rho = 8, \beta = 0.4', ...
    'Proposed Scheme with B_t = B_r = \infty, \rho = 8,\beta = 0.15', ...
    'Proposed Scheme with B_t = B_r = \infty, \rho = 16,\beta = 0.15')
set(gca, 'xtick', SNR);
