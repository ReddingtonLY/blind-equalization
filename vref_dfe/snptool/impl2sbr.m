function [sbr,stp] = impl2sbr(t,impl,gbaud)
    [n0,n1,n2] = size(impl);
    sbr = zeros(n0,n1,n2);
    stp = cumsum(impl);
    ui  = 1e-9 / gbaud;
    for k1 = 1:n1
	for k2 = 1:n2
	    t0 = t(:,k1,k2);
	    t1 = t0 - ui;
	    sbr(:,k1,k2) = stp(:,k1,k2) - interp1(t0,stp(:,k1,k2),t1,'makima',stp(1,k1,k2));
	end
    end
end
