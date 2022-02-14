function s1pout = cats1p(varargin)

    nf = size(varargin{1},1);

    for k = [1:length(varargin)]
	s2pin{k} = zeros(nf,2,2);
	s2pin{k}(:,2,1) = varargin{k};
	s2pin{k}(:,1,2) = varargin{k};
    end

    s2pout = cats2p(s2pin{:});

    s1pout = s2pout(:,2,1);
end

