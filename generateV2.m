function V = generateV2(S, tt, tau, alpha)
% Make covariance matrix. If tau = 0, use alpha instead
if tau == 0
    tau = -(tt(2) - tt(1)) / log(alpha);
end

S  = S(:);
tt = tt(:);
n  = length(tt);

% Upper triangle only (j > i), then mirror
[i, j] = find(triu(true(n), 0));   % indices du triangle supérieur
vals   = S(i) .* S(j) .* exp(-abs(tt(i) - tt(j)) / tau);

V = zeros(n, n);
V(sub2ind([n n], i, j)) = vals;
V = V + triu(V, 1)';   % recopie le triangle sup → inf