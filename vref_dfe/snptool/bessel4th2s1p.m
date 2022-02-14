function s1p = bessel4th2s1p(f,fc)

    w   = f / fc * 2 * pi / 2.978474891650; %/ 2.9785;
    s   = w * j;
    ssq = s .* s;
    den = s.^4 + 10 * s.^3 + 45 * s.^2 + 105 * s + 105;
    s1p = 105 ./ den;

end
