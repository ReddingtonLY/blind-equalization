function t2p = s2p2t2p(s2p)

    nf   = size(s2p,1);
    t2p  = zeros(nf,2,2);
    t2pk = zeros(2,2);

    for k = 1:nf
	s2pk = shiftdim(s2p(k,:,:),1);

	t2pk(1,1) = - det(s2pk) / s2pk(2,1);
	t2pk(1,2) =   s2pk(1,1) / s2pk(2,1);
	t2pk(2,1) = - s2pk(2,2) / s2pk(2,1);
	t2pk(2,2) =   1         / s2pk(2,1);

	t2p(k,:,:) = shiftdim(t2pk, -1);
    end
end
