function snpout = snp2snpc(fin,snpin,varargin)

    % transform SnP to causal SnP by truncating anti-causal response

    [n0,n1,n2] = size(snpin);

    [inp,tnp,dnp] = snp2impl(fin,snpin,varargin{:});	% transform to time domain

    for k1 = 1:n1
	for k2 = 1:n2
	    inp(find(tnp(:,k1,k2) < 0), k1, k2) = 0;	% truncating anti-causal response
	end
    end

    [fout,snpout0] = impl2snp(inp,tnp);			% back to frequency domain

    if (fin(1) == 0)
	snpout = snpout0;				% return with DC gain
    else
	snpout = snpout0(2:end,:,:);			% return without DC gain
    end

end

