
% set path for public C2M channel data

c2mchpath = struct([]);

c2mchdatadir = '/nfshome/redliu/Documents/MATLAB/chdata/public/C2M/';



% public/C2M/JLchan5/

% 31-July-2019, Jane Lim, Cisco
% 100G Host-to-Module Short Channels
%   http://www.ieee802.org/3/ck/public/adhoc/aug14_19/lim_3ck_adhoc_01_073119.pdf
% 2in, 3in, 4in, & 9in host trace
%   http://www.ieee802.org/3/ck/public/tools/c2m/lim_3ck_adhoc_02_073119.zip

% c2mch114 : JLchan5/lim_3ck_adhoc_02_073119/Channel5a_Smaller_Pad_2inch_trace/
% c2mch115 : JLchan5/lim_3ck_adhoc_02_073119/Channel5b_Smaller_Pad_3inch_trace/
% c2mch116 : JLchan5/lim_3ck_adhoc_02_073119/Channel5c_Smaller_Pad_4inch_trace/
% c2mch117 : JLchan5/lim_3ck_adhoc_02_073119/Channel5d_Smaller_Pad_9inch_trace/

intable = [2 3 4 9];

ch = 114; %90;
for k = [1:4]
    c2mchpath{ch}.dir  = sprintf('%sJLchan5/lim_3ck_adhoc_02_073119/Channel5%c_Smaller_Pad_%dinch_trace/', c2mchdatadir, 'a'+k-1, intable(k));
    c2mchpath{ch}.thru = sprintf('Channel5_thru_small_pad_%dinch.s4p', intable(k));
    for xt = [1:3]
	c2mchpath{ch}.next{xt} = sprintf('Channel5_NEXT%d_small_pad_%dinch.s4p', xt, intable(k));
    end
    for xt = [1:5]
	c2mchpath{ch}.fext{xt} = sprintf('Channel5_FEXT%d_small_pad_%dinch.s4p', xt, intable(k));
    end
    c2mchpath{ch}.port = '[1 3 2 4]';
    ch = ch + 1;
end



