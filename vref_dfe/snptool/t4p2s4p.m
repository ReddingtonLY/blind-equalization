function s4p = t4p2s4p(t4p)

    nf   = size(t4p,1);
    s4p  = zeros(nf,4,4);
    s4pk = zeros(4,4);

    for k = 1:nf
	t4pk = shiftdim(t4p(k,:,:),1);
	t2424 = t4pk([2 4],[2 4]);
	t2424_inv = inv(t2424);
	t1324 = t4pk([1 3],[2 4]);

	s4pk([2 4],1) =   t2424_inv * [1; 0];
	s4pk([2 4],2) = - t2424_inv * t4pk([2 4],1);
	s4pk([2 4],3) =   t2424_inv * [0; 1];
	s4pk([2 4],4) = - t2424_inv * t4pk([2 4],3);
	s4pk([1 3],1) =   t1324 * s4pk([2 4],1);
	s4pk([1 3],2) =   t1324 * s4pk([2 4],2) + t4pk([1 3],1);
	s4pk([1 3],3) =   t1324 * s4pk([2 4],3);
	s4pk([1 3],4) =   t1324 * s4pk([2 4],4) + t4pk([1 3],3);

	s4p(k,:,:) = shiftdim(s4pk, -1);
    end
end
