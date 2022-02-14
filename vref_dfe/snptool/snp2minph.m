function minph = snp2minph(fin,snpin)
    if (fin(1) == 0)
	f = fin;
	snp = snpin;
    else
	dcg = snp2dcg(fin, snpin, 0.005, 0.1);
    	f   = [0; fin];
	snp = [shiftdim(dcg,-1); snpin];
    end
    nf = length(f);
    [n0,n1,n2] = size(snpin);
    minph = zeros(n0,n1,n2);
    for k1 = 1:n1
	for k2 = 1:n2
	    s1p       = snp(:,k1,k2);
	    gainlog   = log(abs(s1p));
	    gainsym   = [gainlog; gainlog(end-1:-1:2)];
	    zeroceps  = ifft(gainsym, 'symmetric');
	    minphceps = [zeroceps(1); zeroceps(2:nf-1) * 2; zeroceps(nf); zeros(nf-2,1)];
	    minphsym  = fft(minphceps);
	    if (fin(1) == 0)
		minph(:,k1,k2) = imag(minphsym(1:nf));
	    else
		minph(:,k1,k2) = imag(minphsym(2:nf));
	    end
	end
    end
end
