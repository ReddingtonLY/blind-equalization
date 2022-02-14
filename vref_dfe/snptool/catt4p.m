function t4pout = catt4p(t4pin, varargin)

    nf = size(t4pin,1);

    t4pout = zeros(nf,4,4);

    for kf = 1:nf
	t4pk = shiftdim(t4pin(kf,:,:), 1);
	for k = 1:length(varargin)
	    t4pk = t4pk * shiftdim(varargin{k}(kf,:,:), 1);
	end
	t4pout(kf,:,:) = shiftdim(t4pk, -1);
    end
end
