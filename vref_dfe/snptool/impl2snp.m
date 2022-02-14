function [f,snp] = impl2snp(inp,tnp)

    [n0,n1,n2] = size(inp);

    n = n0 / 2;

    f   = [0:n]' / diff(tnp([1 n+1],1,1)) / 2;
    snp = zeros(n+1,n1,n2);

    for k1 = 1:n1
	for k2 = 1:n2
	    s1psym = fft(inp(:,k1,k2));
	    snp(:,k1,k2) = s1psym(1:n+1) .* sec2snp(f,tnp(1,k1,k2));
	end
    end

end

