function s2pout = cats2p(varargin)

    nf = size(varargin{1},1);

    s2pout = zeros(nf,2,2);

    d21 = zeros(length(varargin),1);

    for kf = 1:nf

	% check zero-transfer blocks

	for k = 1:length(varargin)
	    d21(k) = abs(varargin{k}(kf,2,1));
	end
	zd21 = find(d21 == 0);

	if (length(zd21) == 0)

	    % there is no zero-transfer block, a certain gain between left and right

	    t2pk = eye(2);
	    for k = 1:length(varargin)
		t2pk = t2pk * shiftdim(s2p2t2p(varargin{k}(kf,:,:)), 1);
	    end
	    s2pout(kf,:,:) = t2p2s2p(shiftdim(t2pk, -1));

	else

	    % there is one or more zero-transfer block, no gain between left and right

	    % check S-parameter of LHS of zero-transfer block

	    if (zd21(1) == 1)
		s2p_lhs = varargin{1}(kf,:,:);
	    else
		t2pk = eye(2);
		for k = 1:zd21(1)-1
		    t2pk = t2pk * shiftdim(s2p2t2p(varargin{k}(kf,:,:)), 1);
		end
		s2p_zt1 = shiftdim([0 1;1 0], -1);
		s2p_zt1(1,1,1) = varargin{zd21(1)}(kf,1,1);
		t2pk = t2pk * shiftdim(s2p2t2p(s2p_zt1));
		s2p_lhs = t2p2s2p(shiftdim(t2pk, -1));
	    end
	    s2pout(kf, 1, 1) = s2p_lhs(1, 1, 1);

	    % check S-parameter of RHS of zero-transfer block

	    if (zd21(end) == length(varargin))
		s2p_rhs = varargin{end}(kf,:,:);
	    else
		s2p_zt2 = shiftdim([0 1;1 0], -1);
		s2p_zt2(1,2,2) = varargin{zd21(end)}(kf,2,2);
		t2pk = shiftdim(s2p2t2p(s2p_zt2));
		for k = zd21(end)+1:length(varargin)
		    t2pk = t2pk * shiftdim(s2p2t2p(varargin{k}(kf,:,:)), 1);
		end
		s2p_rhs = t2p2s2p(shiftdim(t2pk, -1));
	    end
	    s2pout(kf, 2, 2) = s2p_rhs(1, 2, 2);

	end

    end
end
