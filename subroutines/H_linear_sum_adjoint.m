function dHdx = H_linear_sum_adjoint(X,w)
% subroutine for the adjoint of the observation operator (analytical
% solution) for "one observation"
% input:
% X : ensmemble in state space (size: [# state variable in inner domain * # of ens member])
% w : the weights for each state variable in inner domain
% output:
% dHdx: adjoint (size: [# state variable in inner domain * # of ens member])
% 2022/02/25

[dim_inner, np] = size(X);
if length(w) ~= dim_inner
    error('length of the weights should be equal to the inner domain size')
end

dHdx = repmat(reshape(w, [dim_inner 1]), [1 np]); 

end