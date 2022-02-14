function snpout = ind2snp(f,np,L)

    r0  = 50;	% reference impedance

    s   = 2j * pi * f;

    zL  = s * L;

    s11 = zL      ./ (2 * r0 + zL);
    s21 = 2  * r0 ./ (2 * r0 + zL);

    nf = length(f);
    snpout = zeros(nf,np,np);
    snpout(:,1,1) = s11;
    if (np >= 2)
	snpout(:,2,2) = s11;
	snpout(:,2,1) = s21;
	snpout(:,1,2) = s21;
	if (np == 4)
	    snpout(:,[3 4],[3 4]) = snpout(:,[1 2],[1 2]);
	end
    end
end
