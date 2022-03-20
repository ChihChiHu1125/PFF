% Particel Flow Filter main script
% 2022/03/02
% by Chih-Chi Hu 
% chihchi@colostate.edu

%% Namelist
% for experiment setup:
DA_run   = 1;                  % DA_run = 1: run the ensemble with DA cycles
noDA_run = 1;                  % noDA_run = 1: run the ensemble without DA cycles
nt       = 200;                % number of integration timestep
t_start  = 1;                  % the (re)start time of the model
warm_nt  = 1000;               % warm up time for the Lorenz model
gen_ens  = 1;                  % generate_ens = 1: (for the first run) / generate_ens = 0 (for other run)
np       = 30;                 % number of particles

% parameters for the L96 model:
dim = 40;                    % dimension for lorenz 96 model
F   = 8;                       % forcing 
dt  = 0.01;                    % time resolution

% settings for generating the ensemble
Q     = 2*eye(dim);            % background error covariance (only for the initial perturbation)
Q_inv = inv(Q);

% settings for DA/PFF 
% /PFF kernel:
alpha    = 1/np;               % the tuning parameter for the covariance of the kernel
% /PFF iteration:
max_pseudo_step     = 150;      % maximum number of iterations
eps_init            = 5e-2;    % initial pseudo-timestep size (learning rate)  
stop_cri            = 1e-3;    % when the norm(grad_KL(s)) < stop_cri, the PFF iteration stops
stop_cri_percentage = 0.05;    % when norm(grad_KL(s)) < stop_cri_percentage*norm(grad_KL(1)), stop the iteration
min_learning_rate   = 1e-5;    % when learning rate < min_learning rate, stop the iteration
% /PFF prior assumption:
io_local       = 1;            % 0: no localization / 1: localization (1 is recommended)
r_influ        = 4;            % localization "radius" (if io_local == 1)
io_gauss_prior = 1;            % 0 = gaussian mixture prior, 1 = gaussain prior (1 is recommended)
inflation_fac  = 1.25;         % inflation factor for prior covariance
tune_C  = 5/inflation_fac;     % make the covariance of the component of gaussian mixture larger (if io_gauss_prior == 0)
% /PFF precondition:
pre_cond  = 1;                 % 0 = posterior covariance, 1 = prior covariance (1 is recommended)
% /SVD
cond_num  = -5;                % condition number for SVD

% settings for observation (note that the observation operator is defined in the subroutine)
da_intv  = 20;                 % obs frequency (how many timesteps between each observation to be assimilated)

% the inner domain for obs (can be manually changed):
% linear identity obs

% This observation is just identity obs, but not observing all grid points
% for example, y1=x2, y2=x4,...
obs_den      = 2;                         % observation density (every obs_den -th grids are observed)
obs_input    = [obs_den:obs_den:dim]';    % a column vector, each row is the inner domain for a obs
[ny_obs, n_inner] = size(obs_input);
inner_domain = mat2cell(obs_input, [ones(1,ny_obs)],[1]);
%}

% linear summation obs (an example if you would like to create your own obs):
%{
% this observation is like y = weights.*X(obs_input)
% for example, y1 = 0.25*(x1+x2+x3+x4), y2 = 0.25*(y3+y4+y5+y6), etc
obs_den      = 4;                         % observation density (every obs_den -th grids are observed)
weights      = [0.25 0.25 0.25 0.25];
obs_input    = [1:4;   3:6;   5:8;   7:10;  9:12;  11:14; 13:16; 15:18;
                17:20; 19:22; 21:24; 23:26; 25:28; 27:30; 29:32; 31:34;
                33:36; 35:38; 37:40; 39 40 1 2];    % a column vector, each row is the inner domain for a obs
[ny_obs, n_inner] = size(obs_input);
inner_domain = mat2cell(obs_input, [ones(1,ny_obs)],[obs_den]);
%}

% generate the observation error first:
rng('default')
obs_err   = 0.3;                                    % observation error standard deviation
R         = obs_err^2 * eye(ny_obs);                % observation error covariance
total_obs = nt/da_intv;                             % total observation times
obs_rnd   = mvnrnd(zeros(ny_obs,1),R,total_obs)';   % Gaussian obs error

% declare some important matrices
if gen_ens == 1 
    prior = zeros(dim,np,total_obs);   % the matrix used to save the prior for every DA cycle
end
obs       = zeros(ny_obs, total_obs);  % the observations (will be determined later)

% an important matrix recording the norm(grad_KL) for all iterations at all
% observation times:
norm_grad_KL = zeros(total_obs,max_pseudo_step); 

%% Integrate the L96 model (warm up):

% integration for the control run (the truth)
Xt = zeros(dim, warm_nt + nt);
Xt (:,1) = F*ones(dim,1);  % initial condition (for steady state)
Xt(dim/5:dim/5:dim) = F+1; % perturbed IC (to generate chaotic behavior)

for t=1:warm_nt + nt -1
    Xt(:,t+1) = L96_RK4(Xt(:,t),dt,F);
end

%% Generate the ensemble
if gen_ens == 1
    ctlmean = Xt(:,warm_nt+1) + mvnrnd(zeros(dim,1),eye(dim),1)';

    % initial condition
    X = zeros(dim,np,nt);
    X (:,:,1) = mvnrnd(ctlmean,Q,np)'; 
end

%% Integration of ensemble & Sequential Data Assimilation

if noDA_run == 1 % run the ensemble without DA
    XnoDA = zeros(dim,np,nt);
    XnoDA(:,:,1) = X (:,:,1);
    for t=1:nt-1
        XnoDA(:,:,t+1) = L96_RK4(XnoDA(:,:,t),dt,F); % integration of the model
    end
end

if DA_run == 1 % run the ensemble with DA
    t=t_start;
while t<nt
disp(['start timestep t=',num2str(t+1)])
% the settings for frequency of DA:
if mod(t+1,da_intv) == 0
    io_obs = 1;
    io_pff = 1;
    disp(['start DA at t=',num2str(t+1)])
else
    io_obs = 0;
    io_pff = 0;
end
    
% Step 1 -- run the model    
X(:,:,t+1) = L96_RK4(X(:,:,t),dt,F); % integration of the model

% Step 2 -- Generate the observation
if io_obs == 1
    obs_time = (t+1)/da_intv; % the "n-th" observation time (should be an integer)
    for i=1:ny_obs
        inner_ind       = inner_domain{i};
        obs(i,obs_time) = H_linear( Xt(inner_ind, warm_nt+t+1) ) + obs_rnd(i,obs_time);
%         obs(i,obs_time) = H_linear_sum( Xt(inner_ind, warm_nt+t+1), weights ) + obs_rnd(i,obs_time);
    end
end

% Step 3 -- Particle Flow Filter (PFF)
if io_pff == 1

% calculate the (inflated) prior error covariance B and C=B/np, where C is
% the covariance for each component in Gaussian mixture prior model
X_tmp  = X(:,:,t+1);                     % the prior ensemble matrix
X_mean = repmat( mean(X_tmp,2), [1 np]); % the mean state of particles (size = dim * np)
C      = inflation_fac*(X_tmp-X_mean)* (X_tmp-X_mean)'/(np-1)/np; % note B is the covariance and here C = B/np

    if io_local == 1 % localization for C (recommended)
        tmp = zeros(dim);
            for i=1:3*r_influ
                tmp = tmp + exp(-i^2/r_influ^2)* ( diag(ones(dim-i,1),i) + diag(ones(dim-i,1),-i) + diag(ones(i,1),-(dim-i)) + diag(ones(i,1),dim-i));
            end
        mask = tmp + diag(ones(dim,1));    
        C = C.*mask;
        C_inv = inv_SVD(C,cond_num, 1); % option = 1, calculate the inverse
    elseif io_local == 0 % without localization for C:
        C_inv = inv_SVD(C,cond_num, 1); % option = 1, calculate the inverse
    end
    
B     = C*np;                               % B: the (localized, if io_local ==1) prior error covariance
B_inv = C_inv / np;                         % B_inv: inverse of B
% the following is not needed, but could be useful for some analysis
% HT    = H_linear(mean(X_tmp,2),obs_den,2)'; % the analytical adjoint (or transpose of linearized H) evaluated at state ensemble mean
% P_inv = B_inv + HT/R*HT';                   % P: the Kalman Filter posterior covariance, P_inv: inverse of P
% P     = inv_SVD(P_inv, cond_num, 1);        % P: the Kalman Filter posterior covariance (calculated by SVD)

% precondition matrix qn ("Q"uasi-"N"ewton method):
if pre_cond == 0     % posterior covariance without inflation factor
     P  = inv_SVD(P_inv, cond_num, 1);
     qn = P/inflation_fac; 
elseif pre_cond == 1 % prior covariance without inflation factor (recommended) 
     qn = B/inflation_fac; 
end

% Start the iteration (in the pseudo-time)================================
last_time_da = t+1; % store the timestep for the DA (for debug)

% parameters to initialize the while loop
s=1;    % pseudo time
ct=0;   % a counter for adpative learning rate
norm_grad_KL(obs_time, 1) = 1e8;         % a random large number to initiate the iteration
eps = eps_init*ones(max_pseudo_step,1);  % pseudo timestep size (or learning rate for the first iteration)      

% the initial condition for the PFF cycles:
pseudo_X        = zeros(dim, np, max_pseudo_step); % dimension of state, # of particles, # of iterations in PFF
pseudo_X(:,:,1) = X_mean + (X_tmp - X_mean);
grad_KL         = zeros(dim,np,max_pseudo_step);   % the particle flow (for all iterations)

% the main PFF iterations:
while norm_grad_KL(obs_time,max(s-1,1)) > stop_cri % the stopping criterion for PFF iteration (could use other...)
    
Hx   = zeros(ny_obs, np);        % the ensemble in obs space 
dHdx = zeros(ny_obs, dim, np);   % the adjoint of obs operator

% evaluate the observation operator and its adjoint:
for i=1:ny_obs
    inner_ind   = inner_domain{i};
    
    % observation operator:
    Hx(i,:)     = H_linear(squeeze(pseudo_X(inner_ind ,:,s)));
%     Hx(i,:)     = H_linear_sum(squeeze(pseudo_X(inner_ind ,:,s)), weights);

    % adjoint of the observation operator
    tmp_dHdx    = H_linear_adjoint(squeeze(pseudo_X(inner_ind ,:,s))); % analytical sol for the adjoint
%     tmp_dHdx    = H_linear_sum_adjoint(squeeze(pseudo_X(inner_ind ,:,s)), weights); % analytical sol for the adjoint 
%     tmp_dHdx    = adjoint_pseudoinverse(squeeze(pseudo_X(inner_ind ,:,s)), Hx(i,:), cond_num); % pseudo inverse for the adjoint
    dHdx(i,inner_ind,:) = tmp_dHdx; % fill in the zeros for the adjoint 
end

% Calculate the gradient of posterior for each particles:
tmp_grad_log_post = zeros(dim,np);

p_obs   = zeros(dim,np); % gradient of the log likelihood
p_bkg   = zeros(dim,np); % gradient of the log prior
numer   = zeros(dim,np); % the numerator of the second term in Pulido and van Leeuwen (2019) eqn (32) [for Gaussian mixture prior]
denom   = zeros(1,np);   % the denominator of the second term in Pulido and van Leeuwen (2019) eqn (32) [for Gaussian mixture prior]
expon   = zeros(1,np);   % the exponent (a scalar) [for Gaussian mixture prior]

for i=1:np % the "i" is to iterate for all particles
    if io_gauss_prior == 0 % gaussian mixture prior
    	for j=1:np % calculate the smallest value in the exponent:
            expon(j) = 0.5*(pseudo_X(:,j,1)-pseudo_X(:,i,s))'*(C_inv/tune_C)*(pseudo_X(:,j,1)-pseudo_X(:,i,s));
        end 
        min_expon = min(expon);
    
        for j=1:np % the "j" is for the summation        
            numer(:,i) = numer(:,i) + (pseudo_X(:,j,1)) * exp( - (expon(j)-min_expon) );
            denom(i)   = denom(i) + exp( - (expon(j)-min_expon) );
        end    
        p_bkg(:,i) = -(C_inv/tune_C)*(pseudo_X(:,i,s)-numer(:,i)/denom(i));     
        
    elseif io_gauss_prior == 1 % gaussian prior (recommended)
        p_bkg(:,i) = -B_inv*(pseudo_X(:,i,s)-mean(X_tmp,2));
    end
    
    HT = dHdx(:,:,i)'; 
    p_obs(:,i) = HT/R*( obs(:,obs_time)- Hx(:,i) );    % assuming Gaussian observation error here
    tmp_grad_log_post(:,i) = p_obs(:,i) + p_bkg(:,i);

end

% Calculate the particle flow:
% Note that this version the kernel is applied onto each dimension
tmp_K          = zeros(dim,np,np); % kernel matrix (will not be stored for each iteration)
tmp_grad_K     = zeros(dim,np,np); % gradient of kernel, repelling force (will not be stored for each iteration)
grad_KL(:,:,s) = zeros(dim,np);    % important to make sure grad_KL starts from zero!
for d=1:dim
    for i=1:np
        for j=1:np
            % calculate the pairwise kernel and its gradient:  
            if j>=i
                tmp_K(d,i,j)      = exp(-0.5*(pseudo_X(d,i,s)-pseudo_X(d,j,s))'*(1/(B(d,d)*alpha))*(pseudo_X(d,i,s)-pseudo_X(d,j,s)));   % kernel
                tmp_grad_K(d,i,j) = -tmp_K(d,i,j)/(B(d,d)*alpha)*(pseudo_X(d,j,s)-pseudo_X(d,i,s)); % the gradient of kernel
            else
                tmp_K(d,i,j)      = tmp_K(d,j,i);
                tmp_grad_K(d,i,j) = -tmp_grad_K(d,j,i);
            end
            grad_KL(d,i,s) = grad_KL(d,i,s) + (tmp_K(d,i,j) * tmp_grad_log_post(d,j) + tmp_grad_K(d,i,j))/np;
        
        end % endfor j
    end % endfor i
end % endfor d

norm_grad_KL(obs_time,s) = sqrt( sum(sum(grad_KL(:,:,s).^2))/(dim*np) ); % norm of grad KL
disp(['iteration  ','s=',num2str(s),' norm =',num2str(norm_grad_KL(obs_time,s)/norm_grad_KL(obs_time,1)*100),' eps = ',num2str(eps(s))])

% "adaptive" learning rate algorithm (might use other algorithms in the future):

if s==1
    stop_cri = stop_cri_percentage*norm_grad_KL(obs_time,1); 
    pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps(s)*qn*grad_KL(:,:,s);
    s  = s+1;
    ct = ct+1;
elseif (eps(s) < min_learning_rate)
    disp('[Note] learning rate too small, break!')
    break
elseif (s>=2) && ( norm_grad_KL(obs_time,s) > 1.02*max(norm_grad_KL(obs_time,s-1:s-1)) ) 
    % when the norm(grad_KL) becomes larger, it could be the learning rate is too large!
    eps(s-1:end) = eps(s)/1.5;    % reduce the learning rate
    s  = s-1;                     % go back one iteration
    ct = 0;                       % reset the counter to zero
    disp(['[Note] The eps has been changed to ',num2str(eps(s)),' ... Redo the PFF based on this eps'])
elseif (ct>=7) && ( norm_grad_KL(obs_time,s) <= 1.02*max(norm_grad_KL(obs_time,s-1:s-1)) )
    % when the norm(grad_KL) becomes smaller in some consecutive
    % iterations, try to increase the learning rate (maybe the learning rate is too small!)
    eps(s:end) = eps(s)*1.5;                                         % increase the learning rate
    pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps(s)*qn*grad_KL(:,:,s);  % use the new learning rate
    s  = s+1;
    ct = 0;   % reset the counter to zero    
else
    pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps(s)*qn*grad_KL(:,:,s);
    s  = s+1;
    ct = ct+1;
end   

if s>=max_pseudo_step
    break
end

end % end while (the main PFF iteration)

s_end = s;
prior(:,:,obs_time) = X(:,:,t+1); % save the prior
X(:,:,t+1) = pseudo_X(:,:,s_end); % save the analysis as the IC for next step

end % end if (the loop for PFF)

t=t+1;
end % end for (the loop for time)

end % end if (the loop for DA_run == 1)
