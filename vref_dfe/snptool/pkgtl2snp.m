function snpout = pkgtl2snp(f,np,zp,Zc)

    % parameters in Table 93A-3
    g0  = 0;
    a1  = 1.734e-3;
    a2  = 1.455e-4;
    tau = 6.141e-3;

    r0  = 50;	% reference impedance

    gamma1 = a1 * (1 + j);					% 93A-10
    gamma2 = a2 * (1 - 2j/pi*log(f/1e9)) + 2j * pi * tau;	% 93A-11
    gamma  = g0 + gamma1 * sqrt(f/1e9) + gamma2 .* f/1e9;	% 93A-9
    gamma(f==0) = g0;

    rho = (Zc - 2 * r0) / (Zc + 2 * r0);			% 93A-12

    s11 = rho * (1 - exp(- gamma * 2 * zp)) ./ (1 - rho^2 * exp(- gamma * 2 * zp));	% 93A-13
    s21 = (1 - rho^2) * exp(- gamma * zp)   ./ (1 - rho^2 * exp(- gamma * 2 * zp));	% 93A-14

    nf = length(f);
    snpout = zeros(nf,np,np);
    if (np == 1)
	snpout(:,1,1) = s21;
    else
	snpout(:,1,1) = s11;
	snpout(:,2,1) = s21;
	snpout(:,1,2) = s21;
	snpout(:,2,2) = s11;
	if (np == 4)
	    snpout(:,[3 4],[3 4]) = snpout(:,[1 2],[1 2]);
	end
    end
end
