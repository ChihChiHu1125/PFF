function dHdx = H_linear_adjoint(X)
% subroutine for the adjoint of the observation operator (analytical
% solution) for "one observation"
% input:
% X : ensmemble in state space (size: [# state variable in inner domain * # of ens member])
% output:
% dHdx: adjoint (size: [# state variable in inner domain * # of ens member])
% 2022/02/25

[dim_inner, np] = size(X);
dHdx = ones(dim_inner,np);

end