function D = DFT_set(N, Q)
% 生成对应的DFT码本集合
% N 为天线数量 ULA/UPA
% Q 量化精度
%%
if nargin == 1
    M = N;
else   
    if length(N) == 2
        M = [2^Q 2^Q];
    else
        M = 2^Q;
    end
end
if length(N) == 1 % DFT 
    n = (0 : N - 1);
    phase = (0 : (M - 1)) * 2 * pi / M;
    D = exp(1i * n.' * phase) / sqrt(N);
    
elseif length(N) == 2 % 2D-DFT 
    phase_v = (0 : (M(1) - 1)) * 2 * pi / M(1); % azimuth
    phase_h = (0 : (M(2) - 1)) * 2 * pi / M(2); % elevation
    
    phase_size = length(phase_v) * length(phase_h);
    phase = zeros(phase_size, 2);  % 存储所有(n1, n2)的组合

    for i = 1 : length(phase_h)
        phase(((i -1) * length(phase_v) + 1) : (i * length(phase_v)), :)...
            = ([phase_v; zeros(1, length(phase_h)) + phase_h(1, i)]).';
    end
    
    D = zeros(N(1) * N(2), length(phase(:, 1)));
    n = 0 : (N(1) - 1);
    m = 0 : (N(2) - 1);
    for i = 1 : length(phase(:, 1))
        e_az = exp(1i * n * phase(i, 1)); 
        e_el = exp(1i * m * phase(i, 2));
        D(:, i) = kron(e_az,e_el).' ./ sqrt(N(1) * N(2));  % User codebook
    end
else
    error('error in data structure: Nr or Nt');
end

