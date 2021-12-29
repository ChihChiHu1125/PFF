% Mapping Particle Filter - L96 Model (linear obs operator)

% 2020/08/08: test the vector kernel with Gaussian mixture prior
% 2020/11/11: check performance

%% Settings for the run
DA_run   = 1;     % DA_run = 1: run the ensemble with DA cycles
noDA_run = 1;     % noDA_run = 1: run the ensemble without DA cycles
nt       = 200;  % number of integration timestep
t_start  = 1;     % the (re)start time of the model
warm_nt  = 1000;  % warm up the Lorenz model
gen_ens  = 1;     % generate_ens = 1: (for the first run) / generate_ens = 0 (for other run)
np       = 20;    % number of particles

% parameters for the L96 model:
dim = 100;
F   = 8;
dt  = 0.01;  % time resolution

% settings for generating the ensemble
Q     = 2*eye(dim);   % background error covariance (only for the initial perturbation)
Q_inv = inv(Q);


% settings for data assimilation
da_intv  = 20;              % how many timesteps between each observation to be assimilated (obs frequency)
io_local = 1;               % 1: localization / 0: no localization
alpha    = 1/np;         % the tuning parameter for the covariance of the kernel
r_influ  = 4;

% settings for MPF
stop_cri        = 1e-3;  % stopping criterion for the iteration
max_pseudo_step = 50;   % maximum number of iterations
eps_init        = 5e-2;   % chosen pseudo-timestep parameter 
reduce_eps_cri  = 1;   % test run before assign the value; can be different for different observations
change_of_eps   = 0;
cond_num        = -5;
io_fix_kernel   = 0;    % 0 = variable kernel, 1 = fixed kernel (0 is default)
io_gauss_prior  = 1;    % 0 = gaussian mixture prior, 1 = gaussain prior (0 is default)
pre_cond        = 1;    % 0 = posterior covariance, 1 = prior covariance 
inflation_fac   = 1.25;    % inflation for the kernel posterior
tune_C          = 5/inflation_fac;    % make the covariance of the component of gaussian mixture larger
crt             = 5;    % cluster analysis criterion


% settings for observation (IMPORTANT: see the end of file for the observation operator!)
% H    = @(x) x;        % observational operator
% dHdx = @(x) eye(dim); % the derivative of observaional operator
obs_den  = 4;           % ?›¸??”å¹¾?‹è?Šæ•¸è§?æ¸?
R        = 0.5*eye(dim/obs_den);  % error covariance of observation (modify this based on the observation operator)
ny_obs   = dim/obs_den;       % the dimension of the observation


% generate the observation error first:
rng('default')
total_obs = nt/da_intv; % total number of observations (in time)
obs_rnd   = mvnrnd(zeros(ny_obs,1),R,total_obs)';
if gen_ens == 1 
    prior     = zeros(dim,np,total_obs); % the matrix to save the prior for every DA cycle
end
n_cluster = zeros(total_obs, 1); % record of how many members in each cluster for each DA cycle


%% Integrate the L96 model (warm up):

% integration for the control run (the truth)
Xt = zeros(dim, warm_nt + nt);
Xt (:,1) = F*ones(dim,1); % initial condition (for steady state)
Xt(dim/5:dim/5:dim) = F+1; 

for t=1:warm_nt + nt -1
    Xt(:,t+1) = L9640_RK4(Xt(:,t),dt,F);
end

%% Generate the ensemble
if gen_ens == 1
    rng('default')
    mvnrnd(zeros(dim,1),eye(dim),1)';
    mvnrnd(zeros(dim,1),eye(dim),1)';
    ctlmean = Xt(:,warm_nt+1) + mvnrnd(zeros(dim,1),eye(dim),1)';

    % initial condition
    X = zeros(dim,np,nt);
    rng('default')
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

if io_fix_kernel == 1
% construct the climatology background error covariacne
tmp = zeros(dim);
for i=1:3*r_influ
    tmp = tmp + exp(-i^2/r_influ^2)* ( diag(ones(dim-i,1),i) + diag(ones(dim-i,1),-i) + diag(ones(i,1),-(dim-i)) + diag(ones(i,1),dim-i));
end
mask = tmp + diag(ones(dim,1));
        
static_cov = zeros(dim);
for t=1:nt
    static_cov = static_cov + cov(XnoDA(:,:,t)')/nt;
end
static_cov     = static_cov.*mask/4;
static_cov_inv = inv_SVD(static_cov,cond_num, 1); % the inverse of the BE matrix
end


norm_grad_KL = zeros(total_obs,max_pseudo_step);

if DA_run == 1 % run the ensemble with DA
    t=t_start;
while t<nt
disp(['start timestep t=',num2str(t+1)])
% the settings for frequency of DA:
    if mod(t+1,da_intv) == 0
        io_obs = 1;
        io_mpf = 1;
        disp(['start DA at t=',num2str(t+1)])
    else
        io_obs = 0;
        io_mpf = 0;
    end
    
% Step 1 -- run the model    
X(:,:,t+1) = L96_RK4(X(:,:,t),dt,F); % integration of the model

% Step 2 -- Generate the observation (NO NEED FOR PERTURBATION!)
if io_obs == 1
    obs_time = (t+1)/da_intv; % ç¬¬å¹¾æ¬¡å?Œå?–è?æ¸¬è?‡æ??
    obs      = H_linear( Xt(:,warm_nt + t+1), obs_den, 1) + obs_rnd(:,obs_time); 
    %obs      = H( Xt(:,warm_nt + t+1), obs_den ) ; 
end

% Step 3 -- Mapping Particle Filter
if change_of_eps == 0
    eps = eps_init;
end

kernel_mean = zeros(max_pseudo_step,1);
spd_up_iter1 = 0;
spd_up_iter2 = 0;
spd_up_iter3 = 0;

if io_mpf == 1
norm_grad_KL(obs_time,:) = 10;    % initial random value (for initiating the iteration)

% To calculate the estimated covariance (C) of each Gaussian mixture and the
% inverse of the estimated covariance (C^{-1}) by SVD: ===================
X_tmp     = X(:,:,t+1);   % the covariance matrix for the prior
X_mean    = repmat( mean(X_tmp,2), [1 np]); % the mean of the state of particles (with size = dim * np)
C         = inflation_fac*(X_tmp-X_mean)* (X_tmp-X_mean)'/(np-1)/np; % the estimated covariance C = B/np

    if io_local == 1 % localization for C:

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

B_inv = C_inv / np;   % the localized background error covariance
HT    = H_linear(mean(X_tmp,2),obs_den,2)';
P_inv = B_inv + HT/R*HT'; % the Kalman Filter posterior covariance
B     = C*np;
%P     = inv_SVD(P_inv, cond_num, 1);

% calculate the inverse of C by SVD:
if pre_cond == 0
     P  = inv_SVD(P_inv, cond_num, 1);
     qn = P/inflation_fac; % (posterior error covariance)
elseif pre_cond == 1
     qn = C*np/inflation_fac; % the quasi-newton method (background error covariance)
end


% Start the iteration (in the pseudo-time)================================
last_time_da = t+1; % store the timestep for the DA (just for debug purpose!)
s=1; % pseudo timestep
eps = eps_init;
grad_KL    = zeros(dim,np,max_pseudo_step);

% the initial condition for the MPF cycles:
pseudo_X        = zeros(dim, np, max_pseudo_step);
pseudo_X(:,:,1) = X_mean + (X_tmp - X_mean);
io_update = ones(dim/obs_den,max_pseudo_step);
jo_update = ones(dim,max_pseudo_step);
ct=0;

while norm_grad_KL(obs_time,max(s-1,1)) > stop_cri

p_obs   = zeros(dim,np); % the first term in eqn (32)
p_bkg   = zeros(dim,np); % the second term in eqn (32)
numer   = zeros(dim,np); % the numerator of the second term in eqn (32)
denom   = zeros(1,np);   % the denominator of the second term in eqn (32)
expon   = zeros(1,np);   % the exponent (a scalar)

% Calculate the gradient of posterior for each particles:
tmp_grad_log_post = zeros(dim,np);


for i=1:np % the "i" is to iterate for all particles
    if io_gauss_prior == 0 % gaussian mixture prior
        
    	for j=1:np % calculate the smallest value in the exponent:
            expon(j) = 0.5*(pseudo_X(:,j,1)-pseudo_X(:,i,s))'*(C_inv/tune_C)*(pseudo_X(:,j,1)-pseudo_X(:,i,s));
        end
    
        min_expon = min(expon);
    
        for j=1:np % the "j" is for the summation        
            numer(:,i) = numer(:,i) + (pseudo_X(:,j,1)) * exp( - (expon(j)-min_expon) );%??†å??
            denom(i)   = denom(i) + exp( - (expon(j)-min_expon) );%??†æ??
        end
        
        p_bkg(:,i) = -(C_inv/tune_C)*(pseudo_X(:,i,s)-numer(:,i)/denom(i));
        
        
    elseif io_gauss_prior == 1 % gaussian prior
        p_bkg(:,i) = -B_inv*(pseudo_X(:,i,s)-mean(X_tmp,2));
    end
    
    HT = H_linear(pseudo_X(:,i,s),obs_den,2)';
    p_obs(:,i) = HT/R*(obs-H_linear(pseudo_X(:,i,s),obs_den,1));
    tmp_grad_log_post(:,i) = p_obs(:,i) + p_bkg(:,i);


end

% Calculate the particle flow:
% Note that this version the kernel is applied onto each dimension
tmp_K          = zeros(dim,np,np);
tmp_grad_K     = zeros(dim,np,np);
grad_KL(:,:,s) = zeros(dim,np);


for d=1:dim
for i=1:np
    for j=1:np
        if io_fix_kernel == 0 % variable kernel
%            tmp_K(i,j)        = 1/(1+(pseudo_X(:,i,s)-pseudo_X(:,j,s))'*(B_inv/alpha)*(pseudo_X(:,i,s)-pseudo_X(:,j,s)));
%            tmp_grad_K(:,i,j) = -2*tmp_K(i,j)^2/alpha*(B_inv)*(pseudo_X(:,j,s)-pseudo_X(:,i,s)); % the gradient of variable kernel
%           tmp_K(d,i,j)        = exp(-0.5*(pseudo_X(d,i,s)-pseudo_X(d,j,s))'*(1/(P(d,d)*alpha))*(pseudo_X(d,i,s)-pseudo_X(d,j,s)));   %variable kernel
%           tmp_grad_K(d,i,j) = -tmp_K(d,i,j)/(P(d,d)*alpha)*(pseudo_X(d,j,s)-pseudo_X(d,i,s)); % the gradient of variable kernel
           if j>=i
               tmp_K(d,i,j)      = exp(-0.5*(pseudo_X(d,i,s)-pseudo_X(d,j,s))'*(1/(B(d,d)*alpha))*(pseudo_X(d,i,s)-pseudo_X(d,j,s)));   %variable kernel
               tmp_grad_K(d,i,j) = -tmp_K(d,i,j)/(B(d,d)*alpha)*(pseudo_X(d,j,s)-pseudo_X(d,i,s)); % the gradient of variable kernel
           else
               tmp_K(d,i,j)      = tmp_K(d,j,i);
               tmp_grad_K(d,i,j) = -tmp_grad_K(d,j,i);
           end

%           tmp_K(d,i,j)        = exp(-0.5*(pseudo_X(d,i,s)-pseudo_X(d,j,s))'*(P_inv(d,d)/alpha)*(pseudo_X(d,i,s)-pseudo_X(d,j,s)));   %variable kernel
%           tmp_grad_K(d,i,j) = -tmp_K(d,i,j)*P_inv(d,d)/(alpha)*(pseudo_X(d,j,s)-pseudo_X(d,i,s)); % the gradient of variable kernel
%           tmp_K(d,i,j)        = exp(-0.5*(pseudo_X(d,i,s)-pseudo_X(d,j,s))'*(1/0.01)*(pseudo_X(d,i,s)-pseudo_X(d,j,s)));   %variable kernel
%           tmp_grad_K(d,i,j) = -tmp_K(d,i,j)/(0.01)*(pseudo_X(d,j,s)-pseudo_X(d,i,s)); % the gradient of variable kernel
        end
        
        grad_KL(d,i,s) = grad_KL(d,i,s) + (tmp_K(d,i,j) * tmp_grad_log_post(d,j) + tmp_grad_K(d,i,j))/np;
%         grad_KL(:,i,s) = grad_KL(:,i,s) + ( tmp_K(i,j) * tmp_grad_log_post(:,j) )/np;
    end
%    grad_KL(d,i,s) = grad_KL(d,i,s)/sum(tmp_K(d,i,:));
end
end

    norm_grad_KL(obs_time,s) = sqrt( sum(sum(grad_KL(:,:,s).^2))/(dim*np) );
    disp('   ')
    disp(['mean grad KL =',num2str( sqrt( sum(sum(grad_KL(:,:,s).^2))/(dim*np)) )])
    %kernel_mean(s) = mean(mean(tmp_K,1),2);
    %disp('   ')
    %disp(['kernel mean =',num2str(mean(mean(tmp_K,1),2))])
    %disp('   ')

% old convergence algorithm:
%{
    if (norm_grad_KL(obs_time,s) > reduce_eps_cri)&&(s==1) % DA cycleå¾????8æ­?
        reduce_eps_cri = norm_grad_KL(obs_time,s)
    elseif norm_grad_KL(obs_time,s) > reduce_eps_cri % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps/2;
        s=max(s-8,1);
        disp(['The eps has been changed to ',num2str(eps),' ... Redo the MPF based on this eps'])
    else
        pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*qn*grad_KL(:,:,s);
        %pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*grad_KL(:,:,s);
        s=s+1;
        disp(['DA iteration step ',num2str(s)])
    end
%}

% new convergence algorithm:

    if (norm_grad_KL(obs_time,s) > reduce_eps_cri)&&(s==1) % DA cycleå¾????8æ­?
        reduce_eps_cri = norm_grad_KL(obs_time,s)
    elseif (eps<1e-5) % DA cycleå¾????8æ­?
%        norm_grad_KL(obs_time,s) = 0;
%        s=s+1
        break
        %{
    elseif (sum(io_update(:,s))==0) && (sum(io_update(:,max(s-1,1)))~=0) && (norm_grad_KL(obs_time,s) <= 1.005*max(norm_grad_KL(obs_time,s-1:s-1))) % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps*100;
        pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*grad_KL(:,:,s).*repmat(jo_update(:,s), [1 np]);
        s=s+1;
        disp(['Increase eps to ',num2str(eps)])
        ct=0;
        %}
    elseif norm_grad_KL(obs_time,s) > 10*reduce_eps_cri % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps*0.3;
        s=max(s-1,1);
        disp(['The eps has been changed to ',num2str(eps),' ... Redo the MPF based on this eps'])
        ct=0;
    elseif norm_grad_KL(obs_time,s) > 2.5*reduce_eps_cri % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps*0.5;
        s=max(s-1,1);
        disp(['The eps has been changed to ',num2str(eps),' ... Redo the MPF based on this eps'])
        ct=0;
        
    elseif (s>=2) && (norm_grad_KL(obs_time,s) > 1.001*max(norm_grad_KL(obs_time,s-1:s-1))) % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps*0.8;
        s=max(s-1,1);
        disp(['The eps has been changed to ',num2str(eps),' ... Redo the MPF based on this eps'])
        ct=0;
        
    elseif (ct>=20) && (norm_grad_KL(obs_time,s) <= 1.001*max(norm_grad_KL(obs_time,s-1:s-1))) % DA cycleå¾????8æ­?
        change_of_eps = 1;
        eps = eps*1.4;
        %pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*qn*grad_KL(:,:,s).*repmat(jo_update(:,s), [1 np]);
        pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*qn*grad_KL(:,:,s);
        s=s+1;
        disp(['Increase eps to ',num2str(eps)])
        ct=0;
        
    else
        %pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*qn*grad_KL(:,:,s).*repmat(jo_update(:,s), [1 np]);
        pseudo_X(:,:,s+1) = pseudo_X(:,:,s) + eps*qn*grad_KL(:,:,s);
        s=s+1;
        ct=ct+1;
        disp(['DA iteration step ',num2str(s)])
        
    end
    %{
    var_yy = var(H(pseudo_X(:,:,s),obs_den)');
    var_y1 = var(H(pseudo_X(:,:,1),obs_den)');
    yy_tmp = var_y1./var_yy;
    io_update(:,s) = (yy_tmp) < (var_y1/0.5 + 1);
    %io_update(:,s) = var_yy > 0.8;
    jo_update(obs_den:obs_den:dim,s) = io_update(:,s);
    sum(io_update(:,s))
    %}
   
    s_end = s-1;
%}    

if s>=max_pseudo_step
    change_of_eps = 0;
    break
end

end % end while


prior(:,:,obs_time) = X(:,:,t+1); % save the prior
X(:,:,t+1) = pseudo_X(:,:,s_end);


end % end if (the loop for MPF)

t=t+1;
end % end for (the loop for time)

end % end if (the loop for DA_run == 1)

%% the observation operator
%{
function X_out = H(X_in,obs_den)
    [dim, np] = size(X_in);
    X_out = X_in(obs_den:obs_den:dim,:);
    
end

function X_out = dHdx(X_in,obs_den)
    [dim, np] = size(X_in);
    X_out = zeros(dim/obs_den,dim);
    for i=1:dim/obs_den
        X_out(i,obs_den*i) = 1;
    end
end
%}