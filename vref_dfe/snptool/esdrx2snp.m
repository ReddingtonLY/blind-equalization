function snpout = esdrx2snp(f,np,R0,C0,R2,C2)

    Rref = 50;	% reference impedance

    s   = 2j * pi * f;
    xf0 = find(f == 0);

    zC0 = 1 ./ (s * C0);
    zC2 = 1 ./ (s * C2);

    zin = 1 ./ (1 / R0 + 1 ./ zC0 + 1 ./ (R2 + zC2));

    zin(xf0) = R0;

    gamma2 = (zin - Rref) ./ (zin + Rref);

    s11 = gamma2;
    s21 = (1 + gamma2) .* zC2 ./ (R2 + zC2);

    s21(xf0) = 1 + gamma2(xf0);

    nf = length(f);
    snpout = zeros(nf,np,np);
    snpout(:,1,1) = s11;
    if (np >= 2)
	snpout(:,2,1) = s21;
	snpout(:,2,2) = 0;
	snpout(:,1,2) = 0;
	if (np == 4)
	    snpout(:,[3 4],[3 4]) = snpout(:,[1 2],[1 2]);
	end
    end
end
