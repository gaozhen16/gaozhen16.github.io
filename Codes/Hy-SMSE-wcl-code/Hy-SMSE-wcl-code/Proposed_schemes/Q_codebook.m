function D = Q_codebook(N, Q, q)
% 新的DFT码本方案，可实现过采样
% N 为天线数量
% Q 为量化比特
% a 为虚拟角度细分成程度
%%
choose = 1; %临时代码
if nargin == 1
    q = N * 8;
    Q = N;
elseif nargin == 2
    q = 2^Q * 8;
    Q = 2^Q;
else
%     disp('ok')
    if q < 0
        q = -q; choose = 0;
    end
    q = 2^Q * q;
    Q = 2^Q;

end

RF = DFT_set(N, log2(q)); % 采样后得到的码本
if length(N) == 2
    N_total = N(1) * N(2);
else
    N_total = N;
end 

% 量化
if choose == 1
    ang = angle(RF) + 0.0001; %计算机计算问题
    n = round(ang / 2 / pi * Q);
    RF = exp(1i * n / Q * 2 * pi) / sqrt(N_total);
    D = (unique(RF','rows'))';
elseif choose == 0
    D = RF;
end
