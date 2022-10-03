function [B,C_tmp] = get_B_C(X,Y,n)


           
        beta = [.1 .3 .1 .3; .7 .9 .7 .9];
        beta2 = [.1 .3 .7 .9; .7 .9 .1 .3];
        beta3 = [.4 .6 .1 .3; .4 .6 .7 .9];
        beta4 = [.1 .3 .4 .6; .7 .9 .4 .6];
    




    B_tmp1 = 1.*(X>=beta(1,1)) & 1.*(X<=beta(1,2)) & 1.*(Y>=beta(1,3)) & 1.*(Y<=beta(1,4));
    B_tmp2 = 1.*(X>=beta(2,1)) & 1.*(X<=beta(2,2)) & 1.*(Y>=beta(2,3)) & 1.*(Y<=beta(2,4));
    B_tmp3 = 1.*(X>=beta2(1,1)) & 1.*(X<=beta2(1,2)) & 1.*(Y>=beta2(1,3)) & 1.*(Y<=beta2(1,4));
    B_tmp4 = 1.*(X>=beta2(2,1)) & 1.*(X<=beta2(2,2)) & 1.*(Y>=beta2(2,3)) & 1.*(Y<=beta2(2,4));
    
    C_tmp1 = 1.*(X>=beta3(1,1)) & 1.*(X<=beta3(1,2)) & 1.*(Y>=beta3(1,3)) & 1.*(Y<=beta3(1,4));
    C_tmp2 = 1.*(X>=beta3(2,1)) & 1.*(X<=beta3(2,2)) & 1.*(Y>=beta3(2,3)) & 1.*(Y<=beta3(2,4));
    C_tmp3 = 1.*(X>=beta4(1,1)) & 1.*(X<=beta4(1,2)) & 1.*(Y>=beta4(1,3)) & 1.*(Y<=beta4(1,4));
    C_tmp4 = 1.*(X>=beta4(2,1)) & 1.*(X<=beta4(2,2)) & 1.*(Y>=beta4(2,3)) & 1.*(Y<=beta4(2,4));
    B_tmp(1,:) = B_tmp1(:);
    B_tmp(2,:) = B_tmp2(:);
    
    B_tmp(3,:) = B_tmp3(:);
    B_tmp(4,:) = B_tmp4(:);
    
    C_tmp(1,:) = C_tmp1(:)/nnz(C_tmp1(:));
    C_tmp(2,:) = C_tmp2(:)/nnz(C_tmp2(:));
    C_tmp(3,:) = C_tmp3(:)/nnz(C_tmp3(:));
    C_tmp(4,:) = C_tmp4(:)/nnz(C_tmp4(:));





B = sparse(double(B_tmp'));

    



end

