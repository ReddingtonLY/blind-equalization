function [z11, u11] = impl2tdr(i11,r0)

    u11 = cumsum(i11);
    z11 = r0 * (1 + u11) ./ (1 - u11);
end
