function snpout = renormsnp(snpin, z0in, z0out)
    nf = size(snpin,1);
    np = size(snpin,2);

    if (length(z0in) == 1)
	z0in = diag(ones(np,1) * z0in);
    elseif ((size(z0in,1) == 1) || (size(z0in,2) == 1))
	z0in = diag(z0in);
    end
    sqrt_z0in = z0in .^ (1/2);

    if (length(z0out) == 1)
	z0out = diag(ones(np,1) * z0out);
    elseif ((size(z0out,1) == 1) || (size(z0out,2) == 1))
	z0out = diag(z0out);
    end
    sqrt_y0out = inv(z0out) .^ (1/2);

    one_np = eye(np);

    snpout = zeros(nf,np,np);
    for k = [1:nf]
	Skin  = shiftdim(snpin(k,:,:), 1);
	Zk    = sqrt_z0in * (one_np + Skin) * inv(one_np - Skin) * sqrt_z0in;
	Skout = (sqrt_y0out * Zk * sqrt_y0out - one_np) * inv(sqrt_y0out * Zk * sqrt_y0out + one_np);
	snpout(k,:,:) = shiftdim(Skout, -1);
    end
end
