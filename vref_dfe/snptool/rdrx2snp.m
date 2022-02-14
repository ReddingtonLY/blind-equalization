function snpout = rdrx2snp(f,np,rd)

    r0  = 50;	% reference impedance

    gamma2 = (rd - r0) / (rd + r0);	% 93A-17

    s11 = gamma2;			% 93A-18
    s21 = 1 + gamma2;
    s12 = 1 - gamma2;

    nf = length(f);
    snpout = zeros(nf,np,np);
    snpout(:,1,1) = s11;
    snpout(:,2,1) = s21;
    snpout(:,1,2) = s12;
    if (np == 4)
	snpout(:,[3 4],[3 4]) = snpout(:,[1 2],[1 2]);
    end
end
