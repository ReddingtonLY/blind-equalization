function snpout = sec2snp(f,secin)
    nf = length(f);
    np = length(secin);
    snpout = zeros(nf,np,np);
    if (np == 1)
	snpout(:,1,1) = exp(- j * 2 * pi * f * secin(1));
    else
	snpout(:,2,1) = exp(- j * 2 * pi * f * secin(1));
	snpout(:,1,2) = exp(- j * 2 * pi * f * secin(2));
	if (np > 2)
	    snpout(:,4,3) = exp(- j * 2 * pi * f * secin(3));
	    snpout(:,3,4) = exp(- j * 2 * pi * f * secin(4));
	end
    end
end
