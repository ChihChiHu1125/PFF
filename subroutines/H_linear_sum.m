function Hx = H_linear_sum(X,w)
% observation operator for linear weighted average obs
% input:
% X : ensmemble in state space (size: [# state variable in inner domain * # of ens member])
% w : the weights for each state variable in inner domain
% output:
% Hx: ensemble in the observation space (size: [# of ens member])
% 2022/02/25

[dim_inner, np] = size(X);    % np: # of ens member
if length(w) ~= dim_inner
    error('length of the weights should be equal to the inner domain size')
end

W  = repmat(reshape(w, [dim_inner 1]), [1 np]); 
Hx = sum(W.*X,1);

end