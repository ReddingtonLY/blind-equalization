function dcg = snp2dcg(f,snpin,prelen,dcintlen)
    prepts = round(2 * f(end) / diff(f(1:2))  * prelen);
    if (abs(f(1)) < 1e3)
	dcg = shiftdim(real(snpin(1,:,:)), 1);
	return
    end
    [n0,n1,n2] = size(snpin);
    ndc = round(n0 * 2 * dcintlen);
    dcg = zeros(n1,n2);
    for k1 = 1:n1
	for k2 = 1:n2
	    s1p  = snpin(:,k1,k2);
	    phs  = unwrap(angle(s1p));
%	    phs  = angle(s1p);
	    dly  = (- phs(end) / pi - prepts) / 2 / f(end);
	    s1p  = s1p .* exp(j * 2 * pi * f * dly);
	    conjsym = zeros(n0*2,1);
	    conjsym(2:n0+1,1) = s1p;
	    conjsym(n0+2:end,1) = conj(s1p(end-1:-1:1));
	    impldbl = ifft(conjsym, 'symmetric');
	    dcg(k1,k2) = - sum(impldbl(end-ndc+1:end)) / dcintlen;
	end
    end
end
