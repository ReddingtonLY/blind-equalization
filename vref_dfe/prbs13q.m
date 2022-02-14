
function [symbols] = prbs13q(pmdlane,swap,precode,graycode)
% default
% swap = 0
% precode = 0 ---- 94.2.2.6 Precoding
% graycode = 1
if nargin == 1
    swap = 0;
    precode = 0;
    graycode = 1;
elseif nargin == 2
    precode = 0;
    graycode = 1;
elseif nargin == 3
    graycode = 1;
end
%---PRBS13q-----%
%PMD lane
num = (2^13-1)*2;
symbols = zeros(1,num/2);
registers = zeros(1,13);
output = zeros(1,num);
% initial set 
% Table 94-11-PRBS13 seeds and initial output
if pmdlane == 0
    registers = [0 0 0 0 0 1 0 1 0 1 0 1 1];
elseif pmdlane == 1
    registers = [0 0 1 1 1 0 1 0 0 0 0 0 1];
elseif pmdlane == 2
    registers = [1 0 0 1 0 0 0 1 0 1 1 0 0];
elseif pmdlane == 3
    registers = [0 1 0 0 0 1 0 0 0 0 0 1 0];
elseif pmdlane == 4 % credo
    registers = [1 1 1 1 1 1 1 1 1 1 1 1 1];
else 
    warning('exceed the maximum lane number!');
end    
    



% generator  polynomial 94-3
% G(X) = 1 + x + x^2 + x^12 + x^13;
% reg1 + reg2 + reg12 + reg13 (output)
for ii=1:1:num
    output(ii)=xor(registers(13),xor(registers(12),xor(registers(1),registers(2))));
    
    for jj=13:-1:2
        registers(jj)=registers(jj-1);
    end
    registers(1)=output(ii);
    
end
binary_op = reshape(output,2,[])';
if graycode && ~swap && ~precode
   %% only graycode 
    for ii = 1:1:num/2
        temp = binary_op(ii,:);
        if temp == [0 0]
            symbols(ii) = -3;
        elseif temp == [0 1]
            symbols(ii) = -1;
        elseif temp == [1 1];
            symbols(ii) = 1;
        else 
            symbols(ii) = 3;
        end
    end        
elseif ~graycode && ~swap && ~precode 
    %% no graycode
    
    for ii = 1:1:num/2
        if binary_op(ii,:) == [0 0]
            symbols(ii) = -3;
        elseif binary_op(ii,:) == [0 1]
            symbols(ii) = -1;
        elseif binary_op(ii,:) == [1 0]
            symbols(ii) = 1;
        else
            symbols(ii) = 3;
        end
    end
    
elseif ~graycode && swap && ~precode 
    %% only swap
    for ii = 1:1:num/2
        if binary_op(ii,:) == [0 0]
            symbols(ii) = -3;
        elseif binary_op(ii,:) == [0 1]
            symbols(ii) = 1;
        elseif binary_op(ii,:) == [1 0]
            symbols(ii) = -1;
        else
            symbols(ii) = 3;
        end
    end
    
elseif graycode && swap && ~precode 
    %% swap then gray code
    for ii = 1:1:num/2
        if binary_op(ii,:) == [0 0]
            symbols(ii) = -3;
        elseif binary_op(ii,:) == [0 1]
            symbols(ii) = 3;
        elseif binary_op(ii,:) == [1 0]
            symbols(ii) = -1;
        else
            symbols(ii) = 1;
        end
    end    
end
    


end
% filename = sprintf('prbs13qlane%d.csv',pmdlane);
% fid = fopen(filename,'wt');
% fprintf(fid,'%d',output);
% fclose(fid)
