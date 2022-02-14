function [s2pa, s2pb] = divs2p(fin, s2pin, atob)

    nf = size(s2p,1);

    s2pa = zeros(nf,2,2);
    s2pb = zeros(nf,2,2);

    [i2p,t2p,d2p,s2p,f] = snp2impl(fin, s2pin);

    Td = (d2p(2,1) + d2p(1,2)) / 2;

    [i11, i12_34_56_77_65_54_21, i11b] = divref(i2p(:,1,1), t2p(:,1,1), Td);
    [i88, i87_65_43_22_34_56_78, i88b] = divref(i2p(:,2,2), t2p(:,2,2), Td);

    i12_34_56_78 = i2p(:,2,1) .* (t2p(:,2,1) <= Td * 2);
    i87_65_43_21 = i2p(:,1,2) .* (t2p(:,1,2) <= Td * 2);

end

function [i11f,i11m,i11b] = divref(i11,t11,Td)

    n  = length(i11);
    nf = find(t11 > Td,   1) - 1;
    nm = find(t11 > Td*3, 1) - 1;

    u11f = cumsum(i11(1:nf));
    u11f = [u11f; u11f(nf) * [n-nf-1:-1:0]'/(n-nf)];
    i11f = [u11f(1); diff(u11f)];

    i11b = i11 - i11f;

    u11m = cumsum(i11b(1:nm));
    u11m = [u11m; u11m(nm) * [n-nm-1:-1:0]'/(n-nm)];
    i11m = [u11m(1); diff(u11m)];

    i11b = i11b - i11m;

end

