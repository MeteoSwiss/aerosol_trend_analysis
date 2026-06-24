function S = generateS2(time, data, gd, window)

N = length(time);
S = zeros(1, N);

% Build index ranges vectorized using cumsum trick
i = (1:N)';
i1 = max(i - window, 1);
i2 = min(i + window, N);

% Compute std for each window
for ix = 1:N
    S(ix) = std(data(i1(ix):i2(ix)));
end

% Fix NaNs
nan_idx = find(isnan(S));
for k = nan_idx
    if k > 1
        S(k) = S(k-1);
    else
        S(k) = std(data(gd));
    end
end