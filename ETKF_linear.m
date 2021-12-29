% ETKF for Lorenz 96
% linear observation
% 2020/5/18
% 2020/11/12 copy

%% Settings for the run
DA_run   = 1;     % DA_run = 1: run the ensemble with DA cycles
noDA_run = 1;     % noDA_run = 1: run the ensemble without DA cycles
nt       = 200;    % number of integration timestep
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
da_intv  = 20;               % how many timesteps between each observation to be assimilated
io_local = 1;               % 1: localization / 0: no localization
r_influ  = 4;
inflation_fac = 1.25;            % alpha: inflation factor

% settings for ETKF
cond_num        = -5;

% settings for observation (IMPORTANT: see the end of file for the observation operator!)
% H    = @(x) x;        % observational operator
% dHdx = @(x) eye(dim); % the derivative of observaional operator
obs_den  = 4;           % 相隔幾個變數觀測
R        = 0.5*eye(dim/obs_den);  % error covariance of observation (modify this based on the observation operator)
R_inv    = inv(R);
ny_obs   = dim/obs_den;       % the dimension of the observation
grid_obs = obs_den:obs_den:dim;           % the grids that are observed

DR_inv_vector = zeros(dim/obs_den, dim); % vector form of the diagonal matrix
if io_local == 1
    for i=1:dim
        % define the correlation matrix (the correlation is done on observation covariance!)
          tmp_1   = exp(-((i-grid_obs)/r_influ).^2);
          tmp_2   = exp(-((i+dim-grid_obs)/r_influ).^2);
          tmp_3   = exp(-((i-dim-grid_obs)/r_influ).^2);
          tmp_12  = max(tmp_1, tmp_2);
          cor_mat = max(tmp_12, tmp_3);
          cor_mat(cor_mat<1e-5) = 0;
%           cor_mat(cor_mat>=1e-5) = 1;
          DR_inv_vector(:,i) = cor_mat.*diag(R_inv)'; % DR_inv: D * R_inv (D is the correlation matrix)
    end
else
    DR_inv_vector = ones(dim/obs_den, dim)*2; % remember to change this if not using constant obs err covariance!
end

% generate the observation error first:
rng('default')
total_obs = nt/da_intv; % total number of observations (in time)
obs_rnd   = mvnrnd(zeros(ny_obs,1),R,total_obs)';
prior     = zeros(dim, np, total_obs);

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

if DA_run == 1 % run the ensemble with DA
    t=t_start;
while t<nt
disp(['start timestep t=',num2str(t+1)])
% the settings for frequency of DA:
    if mod(t+1,da_intv) == 0
        io_obs = 1;
        io_da = 1;
        disp(['start DA at t=',num2str(t+1)])
    else
        io_obs = 0;
        io_da = 0;
    end
    
% Step 1 -- run the model    
X(:,:,t+1) = L96_RK4(X(:,:,t),dt,F); % integration of the model

% Step 2 -- Generate the observation (NO NEED FOR PERTURBATION!)
if io_obs == 1
    obs_time = (t+1)/da_intv; % 第幾次同化觀測資料
    obs      = H( Xt(:,warm_nt + t+1), obs_den ) + obs_rnd(:,obs_time); 
    %obs      = H( Xt(:,warm_nt + t+1), obs_den ) ; 
end

% Step 3 -- Ensemble Transformation Kalman Filter
xa_mean = zeros(dim,1);
xa_pert = zeros(dim,np);

if io_da == 1
% construct the ensemble perturbation matrix:
X_tmp     = X(:,:,t+1);        % the full ensemble matrix (mean+perturbation) 
xf_mean   = mean(X_tmp,2);     % the mean of the perturbation matrix
XX        = sqrt(inflation_fac)*(X_tmp - repmat(xf_mean, [1 np]))/sqrt(np-1); % the prior perturbation matrix

% ensemble observation matrix:
% HH = dHdx(xf_mean, obs_den);
% YY = HH*XX;
tmp_y      = H(sqrt(inflation_fac)*(X_tmp - repmat(xf_mean, [1 np]))+repmat(xf_mean, [1 np]), obs_den);
tmp_y_mean = mean(tmp_y,2); 
YY = (tmp_y - repmat(tmp_y_mean, [1 np]))/sqrt(np-1);

% localization for ETKF (LETKF): the update is done grid-by-grid
for i=1:dim
DR_inv   = diag(DR_inv_vector(:,i));

% the "square" of the transformation matrix (TT' = inv(I+Y'*inv(R)*Y)):
T_sq_inv = eye(np) + YY'*DR_inv*YY;
T_sq     = inv_SVD(T_sq_inv, cond_num, 1); % the inverse of T_square
TT       = inv_SVD(T_sq,     cond_num, 2); % the square root of T_square

% update the mean and the perturbation matrix for the ensemble:
% xa_mean(i) = xf_mean(i) + XX(i,:)*T_sq*YY'*DR_inv*(obs-H(xf_mean,obs_den));
xa_mean(i) = xf_mean(i) + XX(i,:)*T_sq*YY'*DR_inv*(obs-tmp_y_mean);
xa_pert(i,:) = XX(i,:)*TT;

% if mod(i,100)==0
% disp(i)
% end

end

% updated ensemble:
prior(:,:,obs_time) = X(:,:,t+1);
X(:,:,t+1)= repmat(xa_mean, [1 np]) + xa_pert*sqrt(np-1);
% X_new = repmat(xa_mean, [1 np]) + xa_pert;
end % end if (the loop for ETKF)

t=t+1;
end % end for (the loop for time)

end % end if (the loop for DA_run == 1)

%% the observation operator

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