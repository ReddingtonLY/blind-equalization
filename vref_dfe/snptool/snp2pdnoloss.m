function [pdnoloss,pd,phnoloss,ph] = snp2pdnoloss(f,snpin)
    dcg   = snp2dcg(f, snpin, 0.005, 0.1);
    minph = snp2minph(f, snpin);

    if (f(1) == 0)
	f0 = f;
    else
    	f0 = [0; f];
    end
    nf0 = length(f0);
    [n0,n1,n2] = size(snpin);
    ph       = zeros(n0,n1,n2);
    phnoloss = zeros(n0,n1,n2);
    pd       = zeros(n0,n1,n2);
    pdnoloss = zeros(n0,n1,n2);
    for k1 = 1:n1
	for k2 = 1:n2
	    if (dcg(k1,k2) >= 0)
		ph(:,k1,k2)       = unwrap(angle(snpin(:,k1,k2)));
		phnoloss(:,k1,k2) = unwrap(angle(snpin(:,k1,k2) .* exp(- j * minph(:,k1,k2))));
		if (f(1) == 0)
		    pd(2:end,k1,k2)       = - ph(2:end,k1,k2)       ./ f(2:end) / 2 / pi;
		    pdnoloss(2:end,k1,k2) = - phnoloss(2:end,k1,k2) ./ f(2:end) / 2 / pi;
		else
		    pd(1:end,k1,k2)       = - ph(1:end,k1,k2)       ./ f(1:end) / 2 / pi;
		    pdnoloss(1:end,k1,k2) = - phnoloss(1:end,k1,k2) ./ f(1:end) / 2 / pi;
		end
	    else
		ph(:,k1,k2)       = unwrap(angle(- snpin(:,k1,k2))) - pi;
		phnoloss(:,k1,k2) = unwrap(angle(- snpin(:,k1,k2) .* exp(- j * minph(:,k1,k2)))) - pi;
		if (f(1) == 0)
		    pd(2:end,k1,k2)       = - (ph(2:end,k1,k2)       + pi) ./ f(2:end) / 2 / pi;
		    pdnoloss(2:end,k1,k2) = - (phnoloss(2:end,k1,k2) + pi) ./ f(2:end) / 2 / pi;
		else
		    pd(1:end,k1,k2)       = - (ph(1:end,k1,k2)       + pi) ./ f(1:end) / 2 / pi;
		    pdnoloss(1:end,k1,k2) = - (phnoloss(1:end,k1,k2) + pi) ./ f(1:end) / 2 / pi;
		end
	    end

	end
    end
end
