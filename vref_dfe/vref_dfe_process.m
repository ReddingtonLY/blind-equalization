% VREF calculations
%

    
  

    if Ndfe == 0 || n < Ndfe+3
        dk(n) = ynffeo(n);
    else
        switch(dk_hat(n-1))
            case 1
                if dk_hatMode == 1
                    c_error = dfeCoeff(1,1:1:min(Ndfe,n-1))*dk_hat(n-1:-1:max(n-Ndfe,1)).*vpp;
                else
                    c_error = dfeCoeff(1,1:1:min(Ndfe,n-1))*dk(n-1:-1:max(n-Ndfe,1))';
                end
                %vref(1,:) = vref(1,:) + c_error;
                %tgt(1,:) = tgt(1,:) + c_error;
                dk(n) = ynffeo(n) - c_error;
            case 1/3
                if dk_hatMode == 1
                    c_error = dfeCoeff(2,1:1:min(Ndfe,n-1))*dk_hat(n-1:-1:max(n-Ndfe,1)).*vpp;
                else
                    c_error = dfeCoeff(2,1:1:min(Ndfe,n-1))*dk(n-1:-1:max(n-Ndfe,1))';
                end
                %vref(1,:) = vref(1,:) + c_error;
                %tgt(1,:) = tgt(1,:) + c_error;
                dk(n) = ynffeo(n) - c_error;
            case -1/3
                if dk_hatMode == 1
                    c_error = dfeCoeff(3,1:1:min(Ndfe,n-1))*dk_hat(n-1:-1:max(n-Ndfe,1)).*vpp;
                else
                    c_error = dfeCoeff(3,1:1:min(Ndfe,n-1))*dk(n-1:-1:max(n-Ndfe,1))';
                end
                %vref(1,:) = vref(1,:) + c_error;
                %tgt(1,:) = tgt(1,:) + c_error;
                dk(n) = ynffeo(n) - c_error;
            case -1
                if dk_hatMode == 1
                    c_error = dfeCoeff(4,1:1:min(Ndfe,n-1))*dk_hat(n-1:-1:max(n-Ndfe,1)).*vpp;
                else
                    c_error = dfeCoeff(4,1:1:min(Ndfe,n-1))*dk(n-1:-1:max(n-Ndfe,1))';
                end
                %vref(1,:) = vref(1,:) + c_error;
                %tgt(1,:) = tgt(1,:) + c_error;
                dk(n) = ynffeo(n) - c_error;
        end
        %dk(n) = ynffeo(n);
    end

    vpp_dfe = mean(abs(dk))*3/2;



    





















