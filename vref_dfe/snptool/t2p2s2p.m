function s2p = t2p2s2p(t2p)

    nf   = size(t2p,1);
    s2p  = zeros(nf,2,2);
    s2pk = zeros(2,2);

    for k = 1:nf
	t2pk = shiftdim(t2p(k,:,:),1);

	s2pk(1,1) =   t2pk(1,2) / t2pk(2,2);
	s2pk(1,2) =   det(t2pk) / t2pk(2,2);
	s2pk(2,1) =   1         / t2pk(2,2);
	s2pk(2,2) = - t2pk(2,1) / t2pk(2,2);

	s2p(k,:,:) = shiftdim(s2pk, -1);
    end
end
