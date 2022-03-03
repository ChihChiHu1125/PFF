% diagnostic for PFF 
% RMSE and ensemble spread in state space
% please run PFF_xxx.m first to generate the required information
% 2020/11/12 
% 2022/03/02

%% Define input, parameters:
all_ind = [1:1:dim];
obs_ind = [2:2:dim];
non_obs = setdiff(all_ind,obs_ind);
len_obs   = length(obs_ind);
len_nobs  = length(non_obs);

% Calculate the RMSE and the spread of ensemble
title_name = 'RMSE in state space';

XX = X;
XXt = Xt;
XXnoDA = XnoDA;

%% Calculate RMSE and ensemble spread

% RMSE
ens_mean      = squeeze(mean(XX,2));
ens_mean_noDA = squeeze(mean(XXnoDA,2));
RMSE_noDA   = zeros(nt,1);
RMSE_noDA_o = zeros(nt,1);
RMSE_noDA_n = zeros(nt,1);
RMSE   = zeros(nt,1);
RMSE_o = zeros(nt,1);
RMSE_n = zeros(nt,1);
for t=1:nt
    RMSE_noDA(t)   = sqrt(sum( (ens_mean_noDA(:,t)-XXt(:,warm_nt+t)).^2 )/dim);
    RMSE_noDA_o(t) = sqrt(sum( (ens_mean_noDA(obs_ind,t)-XXt(obs_ind,warm_nt+t)).^2 )/len_obs);
    RMSE_noDA_n(t) = sqrt(sum( (ens_mean_noDA(non_obs,t)-XXt(non_obs,warm_nt+t)).^2 )/len_nobs);
    RMSE(t)        = sqrt(sum( (ens_mean(:,t)-XXt(:,warm_nt+t)).^2 )/dim);
    RMSE_o(t)      = sqrt(sum( (ens_mean(obs_ind,t)-XXt(obs_ind,warm_nt+t)).^2 )/len_obs);
    RMSE_n(t)      = sqrt(sum( (ens_mean(non_obs,t)-XXt(non_obs,warm_nt+t)).^2 )/len_nobs);
end

% spread
spread        = zeros(nt,1);
spread_o      = zeros(nt,1);
spread_n      = zeros(nt,1);
spread_noDA   = zeros(nt,1);
spread_noDA_o = zeros(nt,1);
spread_noDA_n = zeros(nt,1);

for t=1:nt
    spread(t)        = mean( std(XX(:,:,t),0,2) );
    spread_o(t)      = mean( std(XX(obs_ind,:,t),0,2) );
    spread_n(t)      = mean( std(XX(non_obs,:,t),0,2) );
    spread_noDA(t)   = mean( std(XXnoDA(:,:,t),0,2) );
    spread_noDA_o(t) = mean( std(XXnoDA(obs_ind,:,t),0,2) );
    spread_noDA_n(t) = mean( std(XXnoDA(non_obs,:,t),0,2) );
end
%}

%% Plot Section

figure;
set(gcf,'color','white')
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

plot(1:nt,RMSE_o,'r',1:nt,RMSE_n,'b','linewidth',3)
hold on
plot(1:nt,RMSE_noDA_o,'r-.',1:nt,RMSE_noDA_n,'b-.','linewidth',3)
legend('RMSE (obs)','RMSE (unobs)', 'RMSE (noDA obs)', 'RMSE (noDA unobs)','fontsize',18,'location','northwest')
%}

set(gca,'fontsize',24)
grid on
axis([0 nt 0 6])
xlabel('timestep')
title(title_name,'interpreter','latex','fontsize',28)
hold off