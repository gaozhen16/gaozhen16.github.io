function [H, H_At, H_Ar, N, H_D, U_at, U_ap] = channel(Nt, Nr, K, ch)
% generate channel model (ULA, UPA) in narrow band: H = H_Ar * H_D * H_At
% channel's basic information: ch{.fc, .Nc, .Np, .max_BS, .max_MS, .sigma, .lambda, .d}
% BS&MS_Path{.elevation, .azimuth, .gain}
% H_At:  AoDs steering vector
% H_Ar:  AoAs steering vector
% N   :  the path of strongest gain 
%% Set parameters
if nargin == 3
    ch.fc = 3e11; % 频率
    ch.Nc = 8;    % 簇
    ch.Np = 10;   % 每个簇中的路径数量
    ch.max_BS = pi; % BS角度范围
    ch.max_MS = pi; % users角度范围
    ch.sigma = 7.5;     % 标准差
    ch.lambda = 3e8 / ch.fc; % length of carrier wave
    ch.d = ch.lambda / 2;    % antenna spacing
end

angle = ch.sigma * pi / 180;%标准差计算出来的间隔角度
k = 2 * pi * ch.d / ch.lambda;
total_Np = ch.Nc * ch.Np;
%% 

switch (length(Nt) + length(Nt))
    case 2  % creat mmWave channel(ULA)
        H = zeros(Nr, Nt, K);
        H_At = zeros(Nt, total_Np, K);
        H_Ar = zeros(Nr, total_Np, K);
        H_D = zeros(total_Np, total_Np, K);
        for i = 1 : K
            N_t = (0 : (Nt - 1))';    % transmiter
            phi = (ch.max_BS - angle) * rand(ch.Nc, 1) * ones(1, ch.Np) + ...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_BS / 2 + angle / 2;
            phi = k * sin(reshape(phi',[1, total_Np]));
            a_BS = N_t * phi;
            At = exp(1i * a_BS) / sqrt(Nt);

            N_r = (0 : (Nr - 1))';    % receiver
            theta = (ch.max_MS - angle) * rand(ch.Nc, 1) * ones(1,ch. Np) + ...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_MS / 2 + angle / 2;
 
            theta = k * (reshape(theta',[1, total_Np]));
            a_MS = N_r * theta;
            Ar = exp(1i * a_MS) / sqrt(Nr);

            D_am = (randn(1, total_Np) + 1i * randn(1, total_Np)) / sqrt(2); % Path gain
            [~, N] = max(D_am);
            D = diag(D_am) * sqrt(Nt * Nr / total_Np);
            h = Ar * D * At';

            H(:, :, i) = h;
            H_At(:, :, i) = At;
            H_Ar(:, :, i) = Ar;
            H_D(:, :, i) = D;
        end
    case 4  % creat mmWave channel(UPA)
        total_Nt = Nt(1) * Nt(2); total_Nr = Nr(1) * Nr(2);
        H = zeros(total_Nr, total_Nt, K);
        H_At = zeros(total_Nt, ch.Nc * ch.Np, K);
        H_Ar = zeros(total_Nr, ch.Nc * ch.Np, K);
        H_D = zeros(ch.Nc * ch.Np, ch.Nc * ch.Np, K);
        U_at = zeros(1, ch.Nc * ch.Np, K);
        U_ap = zeros(1, ch.Nc * ch.Np, K);
        for i = 1 : K
            n_t = (0 : Nt(1) - 1)';     % transmitter
            m_t = (0 : Nt(2) - 1)';
            phi1 = (ch.max_BS - angle) * rand(ch.Nc, 1) * ones(1, ch.Np) +...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_BS / 2 + angle / 2; % azimuth

            phi1 = reshape(phi1', [total_Np, 1]);
            theta1 = ch.max_BS * rand(ch.Nc, 1) * ones(1, ch.Np) +...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_BS / 2 + angle / 2; % elevation
            
            theta1 = reshape(theta1', [total_Np, 1]);
            U_at(:, :, i) = theta1.';
            U_ap(:, :, i) = phi1.';   %->提取角度
            A_t = zeros(total_Nt, total_Np);
            for path = 1 : total_Np
                e_a1 = exp(-1i * k * sin(phi1(path, 1)) * cos(theta1(path, 1)) * n_t); % channel model 2
                e_e1 = exp(-1i * k * sin(theta1(path, 1)) * m_t);
                A_t(:, path) = kron(e_a1, e_e1) / sqrt(total_Nt);
            end
            
            n_r = (0:(Nr(1) - 1))';     % receiver
            m_r = (0:(Nr(2) - 1))';
            phi2 = (ch.max_MS - angle) * rand(ch.Nc, 1) * ones(1, ch.Np) +...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_MS / 2 + angle / 2; % azimuth
            
            phi2 = reshape(phi2', [total_Np, 1]);
            theta2 = (ch.max_MS - angle) * rand(ch.Nc, 1) * ones(1, ch.Np) +...
                angle * (rand(ch.Nc, ch.Np) - 0.5 * ones(ch.Nc, ch.Np)) - ch.max_MS / 2 + angle / 2; % elevation
            
            theta2 = reshape(theta2', [total_Np, 1]);
            A_r = zeros(total_Nr, total_Np);
            for path = 1 : total_Np
                e_a2 = exp(-1i * k * sin(phi2(path, 1)) * cos(theta2(path, 1)) * n_r); % channel model 2
                e_e2 = exp(-1i * k * sin(theta2(path, 1)) * m_r);
                A_r(:,path) = kron(e_a2,e_e2) / sqrt(total_Nr);
            end
            
%             D_am = (ones(1, total_Np) + 1i * ones(1, total_Np)) / sqrt(2); % Path gain
            D_am = (randn(1, total_Np) + 1i * randn(1, total_Np)) / sqrt(2); % Path gain
            [~, N] = max(D_am);
            D = diag(D_am) * sqrt(total_Nt * total_Nr / total_Np);
            
            h = A_r * D * A_t';

            H(:, :, i) = h;
            H_At(:, :, i) = A_t;
            H_Ar(:, :, i) = A_r;
            H_D(:, :, i) = D;
        end
    otherwise, error('error in data structure: Nr or Nt');
end
