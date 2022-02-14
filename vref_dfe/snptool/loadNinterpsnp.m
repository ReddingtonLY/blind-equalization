function [snpout,dnp] = loadNinterpsnp(fname, port_order, freqout, varargin)

    % load SnP file

    [freqin, snpin] = loadsnp(fname, port_order);

    % interpolate at required frequencies

    [snpout, dnp] = interpsnp(freqin, snpin, freqout, varargin{:});

end
