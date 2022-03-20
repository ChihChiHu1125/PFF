function dHdx = adjoint_pseudoinverse(X,Hx,cond_num)
% subroutine for estimating the adjoint of the observation operator (by assuming linear H)
% input:
% X : ensmemble in state space (size: [# state variable in inner domain * # of ens member])
% Hx: ensemble in the observation space (size: [# of ens member])
% cond_num: condition number used for SVD 
% (neglect the singular values that are smaller than 10^(cond_num)*largest singular value)
% output:
% dHdx: adjoint (size: [# state variable in inner domain * # of ens member])
% 2022/02/25 

[dim_inner, np] = size(X);
pseudo_inv_X = inv_SVD(X',cond_num,1);
% pseudo_inv_X = pinv(X',1e-5);
HT           = pseudo_inv_X*reshape(Hx, [np 1]);

dHdx = repmat(HT, [1 np]); 
    
end