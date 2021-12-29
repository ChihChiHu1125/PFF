% diag MPF_v3.4 (diagnose in the observational space)
% 2020/06/09

dim       = 1000;
nt        = 200;
warm_nt   = 1000;
obs_den   = 4;
da_intv   = 20;
total_obs = nt/da_intv;

rng('default')
obs_rnd   = mvnrnd(zeros(ny_obs,1),R,total_obs)';

all_ind = [1:1:dim];
obs_ind = [4:4:dim];
non_obs = setdiff(all_ind,obs_ind);
len_obs   = length(obs_ind);
len_nobs  = length(non_obs);


% Calculate the RMSE and the spread of ensemble
XX =X_sq_fix_kernel_width_5_alpha_20;
XXnoDA = XnoDA;

%% observed variables:

name1 = 'RMSE obs var (obs space) PFF $$ H({\bf{x}}) = x^2_{(4a)} $$';
%name1 = 'RMSE obs var (obs space) PFF $$ H({\bf{x}}) = |x_{(4a)}| $$';
%name1 = 'RMSE obs var (obs space) PFF $$ H({\bf{x}}) = e^{(\frac{x_{(4a)}}{6})} $$';

% exp(x/6)
% YY     = squeeze(exp(XX(obs_ind, :,:)/6));
% YYnoDA = squeeze(exp(XXnoDA(obs_ind,:,:)/6));
% YXt    = exp(Xt(obs_ind,warm_nt+1:end)/6);

% square
YY     = squeeze(XX(obs_ind, :,:)).^2;
YYnoDA = squeeze(XXnoDA(obs_ind,:,:)).^2;
YXt    = Xt(obs_ind,warm_nt+1:end).^2;

% absolute value
% YY     = abs(squeeze(XX(obs_ind,:,:)));
% YYnoDA = abs(XXnoDA(obs_ind,:,:));
% YXt    = abs(Xt(obs_ind,warm_nt+1:end));

% reproduce the observations:
observation = zeros(dim/obs_den, total_obs);
obs_time_ind = [da_intv:da_intv:nt];

% ???
%{
for t=1:total_obs
    observation(:,t) = Xt(obs_ind,warm_nt+da_intv*t).^2 + obs_rnd(:,t); % square
end

err = sqrt(0.5)*ones(total_obs,1); % the standard deviation of the observation error

%{
figure;
set(gcf,'color','white') % è¨­å?•ç‚º?™½?‰²
set(gcf,'units','centimeters','position',[2 2 60 12.5])
plot(squeeze(YY(var_id/obs_den,:,:))','color',[.7 .7 .7],'linewidth',2.5)
hold on
plot(1:1:nt,YXt(var_id/obs_den,:))
errorbar(obs_time_ind, observation(var_id/obs_den,:), err,'s','markersize',8,'markeredgecolor','red','markerfacecolor','red')
hold off
axis([0 1500 0 120])
%}
%}


% RMSE
ens_mean      = squeeze(mean(YY,2));
ens_mean_noDA = squeeze(mean(YYnoDA,2));
RMSE_noDA     = zeros(nt,1);
RMSE          = zeros(nt,1);

for t=1:nt
    RMSE_noDA(t)   = sqrt(sum( (ens_mean_noDA(:,t)-YXt(:,t)).^2 )/len_obs);
    RMSE(t)        = sqrt(sum( (ens_mean(:,t)-YXt(:,t)).^2 )/len_obs);
end

% spread
variance      = zeros(len_obs,nt);
variance_noDA = zeros(len_obs,nt);
spread    = zeros(nt,1);
spread_noDA = zeros(nt,1);

% calculate the variance first:
for t=1:nt
    for d=1:len_obs
        variance(d,t) = var(squeeze(YY(d,:,t)));
        variance_noDA(d,t) = var(squeeze(YYnoDA(d,:,t)));
    end
end

for t=1:nt
    spread(t) = sqrt(sum( variance(:,t) )/len_obs);
    spread_noDA(t) = sqrt(sum( variance_noDA(:,t) )/len_obs);
end

%}

%% unobserved variables:
%{
name1 = 'RMSE unobs var (obs space) LETKF $$ H({\bf{x}}) = x^2_{(4a)} $$';
%name1 = 'RMSE unobs var (obs space) LETKF $$ H({\bf{x}}) = |x_{(4a)}| $$';
%name1 = 'RMSE unobs var (obs space) LETKF $$ H({\bf{x}}) = e^{(\frac{x_{(4a)}}{6})} $$';

% exp(x/6)
%YY     = squeeze(exp(XX(non_obs, :,:)/6));
%YYnoDA = squeeze(exp(XXnoDA(non_obs,:,:)/6));
%YXt    = exp(Xt(non_obs,warm_nt+1:end)/6);

% square
% YY     = squeeze(XX(non_obs, :,:)).^2;
% YYnoDA = squeeze(XXnoDA(non_obs,:,:)).^2;
% YXt    = Xt(non_obs,warm_nt+1:end).^2;

% absolute value
YY     = abs(squeeze(XX(non_obs,:,:)));
YYnoDA = abs(XXnoDA(non_obs,:,:));
YXt    = abs(Xt(non_obs,warm_nt+1:end));



% RMSE
ens_mean      = squeeze(mean(YY,2));
ens_mean_noDA = squeeze(mean(YYnoDA,2));
RMSE_noDA     = zeros(nt,1);
RMSE          = zeros(nt,1);

for t=1:nt
    RMSE_noDA(t)   = sqrt(sum( (ens_mean_noDA(:,t)-YXt(:,t)).^2 )/len_nobs);
    RMSE(t)        = sqrt(sum( (ens_mean(:,t)-YXt(:,t)).^2 )/len_nobs);
end

% spread
variance      = zeros(len_nobs,nt);
variance_noDA = zeros(len_nobs,nt);
spread    = zeros(nt,1);
spread_noDA = zeros(nt,1);

% calculate the variance first:
for t=1:nt
    for d=1:len_nobs
        variance(d,t) = var(squeeze(YY(d,:,t)));
        variance_noDA(d,t) = var(squeeze(YYnoDA(d,:,t)));
    end
end

for t=1:nt
    spread(t) = sqrt(sum( variance(:,t) )/len_nobs);
    spread_noDA(t) = sqrt(sum( variance_noDA(:,t) )/len_nobs);
end

%% ========================================================================
%}

figure;
set(gcf,'color','white') % è¨­å?•ç‚º?™½?‰²
set(gcf,'units','centimeters','position',[2 2 26 17])

% compare of spread vs RMSE
%{
plot(1:nt,spread_noDA,'r','linewidth',3)
hold on
plot(1:nt,RMSE_noDA,'r-.','linewidth',3.5)

plot(1:nt,spread,'b','linewidth',3.5)
plot(1:nt,RMSE,'b-.','linewidth',3)
hold off
legend('spread (noDA)','RMSE (noDA)','spread (DA)','RMSE (DA)')
%}

% compare the RMSE of different subsets:

%plot(1:nt,S_exp_n_r6(1:nt),'color',[.8,.2,.2],'linewidth',2.5)
plot(1:nt,spread,'color',[.8, .2, .2],'linewidth',2.5)
hold on
%plot(1:nt,R_exp_n_r6(1:nt),'color',[.3,.3,.3],'linewidth',2.5)
plot(1:nt,RMSE,'color',[.3 .3 .3],'linewidth',2.5)
plot(1:nt,spread_noDA,'-.','color',[.8, .2, .2],'linewidth',2.5)
plot(1:nt,RMSE_noDA,'-.','color',[.3, .3, .3],'linewidth',2.5)
legend('spread (DA)', 'RMSE (DA)','spread (noDA)', 'RMSE (noDA)','fontsize',13,'location','northwest')
% legend('spread','spread (noDA)','RMSE','RMSE (noDA)','location','northwest')


% Plot the RMSE at the observation times:
%{
plot(obs_time_ind,spread(obs_time_ind),'color',[.8 .2 .2],'linewidth',2.5)
hold on
plot(obs_time_ind,RMSE(obs_time_ind),'color',[.3, .3, .3],'linewidth',2.5)
plot(obs_time_ind,spread_noDA(obs_time_ind),'-.','color',[.8 .2 .2],'linewidth',2.5)
plot(obs_time_ind,RMSE_noDA(obs_time_ind),'-.','color',[.3, .3, .3],'linewidth',2.5)
legend('RMSE (DA)','RMSE (noDA)','fontsize',13,'location','northwest')
hold off
%}

set(gca,'fontsize',24)
grid on

% square
axis([0 nt 0 100])

% exp
% axis([0 nt 0 10])

% abs or linear
% axis([0 nt 0 4])

xlabel('timestep')
title('LETKF real6 unobs (obs space)')
%title(name1,'interpreter','latex','fontsize',28)
%title('ETKF (obs space) $$ H(x) = x^2_{4i} $$ (obs freq = 5)','interpreter','latex','fontsize',24)
%title('MPF $$ H(x) = x^2 $$ obs freq = 20, G prior, variabe kernel','interpreter','latex','fontsize',22)
%title('MPF $$ H(x) = x_{4i-1}x_{4i} $$ obs freq = 20, G prior, fixed kernel','interpreter','latex','fontsize',22)
hold off
%}