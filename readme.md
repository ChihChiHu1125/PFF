# Particle Flow Filter for Lorenz 96 

## MATLAB files:
- *PFF.m*: the main script for running data assimilation experiment
- */diagnostics/*: the output diagnostic file for experiment results (from running PFF.m)
- */subroutines/*: contains Lorenz model, SVD, observation operators and their adjoint

## Tips for changing code (for testing different experiment)

### Change the observation
In order to change the observation, we need to do the following:

First need to tell the MATLAB code what is the **inner domain** for each observation.

The **inner domain** for an observation is defined as *the subset of state variables that are used to evaluate the observation operator.* 
For example: In an Lorenz 96 model with 40 variables (denoted as x1,x2,...,x40), let's assume we have an observation operator:
```
H1(x1,x2,...,x40) = x1+x2
```
Note that H1 can be written as H1(x1,x2), since it only requires information from x1 and x2 (not from x3,x4,...,x40). 
So the inner domain of H1 is {x1,x2}

The following is the default observations in PFF.m:
now consider we would like to observe observe every other grid point (i.e., x2, x4, ..., x40).
In this case, each observation is "identity mapping", which is defined by the subroutine "H_linear" (as below), which basically
takes whatever input is as the output. In this case, the inner domain for the first observation is {x2} , {x4} for the second observation, ... 
{x40} for the 20-th observation. There are 40/2 observations. We define the inner domain by the following lines:

```
dim = 40;
% This observation is just identity obs, but not observing all grid points
% for example, y1=x2, y2=x4,...
obs_den      = 2;                         % observation density (every obs_den -th grids are observed)
obs_input    = [obs_den:obs_den:dim]';    % a column vector, each row is the inner domain for a obs
[ny_obs, n_inner] = size(obs_input);
inner_domain = mat2cell(obs_input, [ones(1,ny_obs)],[1]);
```

Note that obs_input defines the inner domain (specifically, each row is the inner domain for each observation).
Note that we turn obs_input from a matrix to a "cell" inner_domain. By doing this, we allow variable size of inner domain for each observation.

```
function Hx = H_linear(X)
[dim_inner, np] = size(X);    % np: # of ens member
Hx = X;
end
```

We also need to modify the lines where the code evaluate the observation operator, 

```
% Step 2 -- Generate the observation
if io_obs == 1
    obs_time = (t+1)/da_intv; % the "n-th" observation time (should be an integer)
    for i=1:ny_obs
        inner_ind       = inner_domain{i};
        obs(i,obs_time) = H_linear( Xt(inner_ind, warm_nt+t+1) ) + obs_rnd(i,obs_time);
    end
end
```

Just need to make sure we are referring to the right subroutine. For example, here we use "H_linear".
Note that the above lines is how we "generate synthetic observation" from the truth (assuming perfect observation operator).
We also need to evaluate the model equivalence (i.e., the ensemble in the observation space), by the following:

```
% evaluate the observation operator and its adjoint:
for i=1:ny_obs
    inner_ind   = inner_domain{i};
    
    % observation operator:
    Hx(i,:)     = H_linear(squeeze(pseudo_X(inner_ind ,:,s)));

    % adjoint of the observation operator
    tmp_dHdx    = H_linear_adjoint(squeeze(pseudo_X(inner_ind ,:,s))); % analytical sol for the adjoint
    dHdx(i,inner_ind,:) = tmp_dHdx; % fill in the zeros for the adjoint 
end
```

Note here we need the adjoint of the observation operator, which is defined by the subroutine H_linear_adjoint.

If you would like to create your own observations, please remember to both code the "observation operator" (e.g., H_linear)
and its associated "adjoint" (e.g., H_linear_adjoint)

