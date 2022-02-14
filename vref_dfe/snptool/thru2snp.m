function snpout = thru2snp(f,np)

    nf = length(f);
    snpout = zeros(nf,np,np);
    if (np == 1)
	snpout(:,1,1) = 1;
    else
	snpout(:,2,1) = 1;
	snpout(:,1,2) = 1;
	if (np > 2)
	    snpout(:,4,3) = 1;
	    snpout(:,3,4) = 1;
	end
    end
end
