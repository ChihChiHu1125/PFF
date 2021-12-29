function X_out = H_linear(X_in,obs_den,option)
% 2020/11/11: observation operator subroutine
% option = 1; X_out = H(X_in)
% option = 2; X_out = dHdx(X_in)

    [dim, np] = size(X_in);
    if option==1
        X_out = X_in(obs_den:obs_den:dim,:);
    elseif option==2
        X_out = zeros(dim/obs_den,dim);
        
        for i=1:dim/obs_den
            X_out(i,obs_den*i) = 1;
        end
    end
        
    
end
