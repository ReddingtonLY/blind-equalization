function zout = zcap(f,c)

    s    = 2j * pi * f;

    zout = 1 ./ (s * c);

end
