function snpout = zs2snp(f,np,zs)

    r0  = 50;	% reference impedance

    s11 =       zs ./ (zs + 2 * r0);
    s21 =   2 * r0 ./ (zs + 2 * r0);

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
