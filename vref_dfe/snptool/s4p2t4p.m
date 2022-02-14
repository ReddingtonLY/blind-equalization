function t4p = s4p2t4p(s4p)

    nf   = size(s4p,1);
    t4p  = zeros(nf,4,4);
    t4pk = zeros(4,4);

    for k = 1:nf
	s4pk = shiftdim(s4p(k,:,:),1);
	s2413 = s4pk([2 4],[1 3]);
	s2413_inv = inv(s2413);
	s1313 = s4pk([1 3],[1 3]);

	t4pk([2 4],1) = - s2413_inv * s4pk([2 4],2);
	t4pk([2 4],2) =   s2413_inv * [1; 0];
	t4pk([2 4],3) = - s2413_inv * s4pk([2 4],4);
	t4pk([2 4],4) =   s2413_inv * [0; 1];
	t4pk([1 3],1) =   s1313 * t4pk([2 4],1) + s4pk([1 3],2);
	t4pk([1 3],2) =   s1313 * t4pk([2 4],2);
	t4pk([1 3],3) =   s1313 * t4pk([2 4],3) + s4pk([1 3],4);
	t4pk([1 3],4) =   s1313 * t4pk([2 4],4);

	t4p(k,:,:) = shiftdim(t4pk, -1);
    end
end
