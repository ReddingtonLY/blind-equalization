function s4pout = cats4p(varargin)

    nf = size(varargin{1},1);

    s4pout = zeros(nf,4,4);

    d2413 = zeros(length(varargin),1);

    for kf = 1:nf

	% check zero-transfer blocks

	for k = 1:length(varargin)
	    d2413(k) = det(shiftdim(varargin{k}(kf,[2 4],[1 3]), 1));
	end
	zd2413 = find(d2413 == 0);

	if (length(zd2413) == 0)

	    % there is no zero-transfer block, a certain gain between left and right

	    t4pk = eye(4);
	    for k = 1:length(varargin)
		t4pk = t4pk * shiftdim(s4p2t4p(varargin{k}(kf,:,:)), 1);
	    end
	    s4pout(kf,:,:) = t4p2s4p(shiftdim(t4pk, -1));

	else

	    % there is one or more zero-transfer block, no gain between left and right

	    % check S-parameter of LHS of zero-transfer block

	    if (zd2413(1) == 1)
		s4p_lhs = varargin{1}(kf,:,:);
	    else
		t4pk = eye(4);
		for k = 1:zd2413(1)-1
		    t4pk = t4pk * shiftdim(s4p2t4p(varargin{k}(kf,:,:)), 1);
		end
		s4p_zt1 = shiftdim([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0], -1);
		s4p_zt1(1,[1 3],[1 3]) = varargin{zd2413(1)}(kf,[1 3],[1 3]);
		t4pk = t4pk * shiftdim(s4p2t4p(s4p_zt1));
		s4p_lhs = t4p2s4p(shiftdim(t4pk, -1));
	    end
	    s4pout(kf, [1 3], [1 3]) = s4p_lhs(1, [1 3], [1 3]);

	    % check S-parameter of RHS of zero-transfer block

	    if (zd2413(end) == length(varargin))
		s4p_rhs = varargin{end}(kf,:,:);
	    else
		s4p_zt2 = shiftdim([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0], -1);
		s4p_zt2(1,[2 4],[2 4]) = varargin{zd2413(end)}(kf,[2 4],[2 4]);
		t4pk = shiftdim(s4p2t4p(s4p_zt2));
		for k = zd2413(end)+1:length(varargin)
		    t4pk = t4pk * shiftdim(s4p2t4p(varargin{k}(kf,:,:)), 1);
		end
		s4p_rhs = t4p2s4p(shiftdim(t4pk, -1));
	    end
	    s4pout(kf, [2 4], [2 4]) = s4p_rhs(1, [2 4], [2 4]);

	end

    end
end
