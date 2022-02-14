function s1p = buttern2s1p(f,fc,n)

    w   = f / fc;
    s   = w * j;
    ssq = s .* s;
    if (mod(n,2) == 0)
	den = ones(size(f));
    else
	den = s + 1;
    end
    for k = 1:floor(n/2)
	den = den .* (ssq - 2 * cos((2 * k + n - 1) / 2 / n * pi) * s + 1);
    end
    s1p = 1 ./ den;

end
