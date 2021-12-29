% diag MPF_v3.4 (diagnose in the model space)
% 2020/11/12 : this script is used for fast diag of the DA performance (in model sapce)
% RMSE time series of observed and unobserved variables

all_ind = [1:1:dim];
obs_ind = [4:4:dim];
non_obs = setdiff(all_ind,obs_ind);
len_obs   = length(obs_ind);
len_nobs  = length(non_obs);


% Calculate the RMSE and the spread of ensemble
XX = X;
XXt = Xt;
XXnoDA = XnoDA;

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
variance      = zeros(dim,nt);
variance_noDA = zeros(dim,nt);
spread    = zeros(nt,1);
spread_o  = zeros(nt,1);
%spread_n  = zeros(nt,1);
spread_noDA = zeros(nt,1);
spread_noDA_o = zeros(nt,1);
%spread_noDA_n = zeros(nt,1);

% calculate the variance first:
%{
for t=1:nt
    for d=1:dim
        variance(d,t) = var(squeeze(XX(d,:,t)));
        variance_noDA(d,t) = var(squeeze(XXnoDA(d,:,t)));
    end
end

for t=1:nt
    spread(t) = sqrt(sum( variance(:,t) )/dim);
    spread_o(t) = sqrt(sum( variance(obs_ind ,t) )/len_obs);
%    spread_n(t) = sqrt(sum( variance(non_obs ,t) )/len_nobs);
    spread_noDA(t) = sqrt(sum( variance_noDA(:,t) )/dim);
    spread_noDA_o(t) = sqrt(sum( variance_noDA(obs_ind ,t) )/len_obs);
%    spread_noDA_n(t) = sqrt(sum( variance_noDA(non_obs ,t) )/len_nobs);
end

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

%plot(1:nt,RMSE,'k',1:nt,RMSE_o,'r',1:nt,RMSE_n,'b','linewidth',3)
plot(1:nt,RMSE_o,'r',1:nt,RMSE_n,'b','linewidth',3)
hold on
plot(1:nt,RMSE_noDA_o,'r-.',1:nt,RMSE_noDA_n,'b-.','linewidth',3)
legend('RMSE (obs)','RMSE (unobs)', 'RMSE (noDA obs)', 'RMSE (noDA unobs)','fontsize',18,'location','northwest')

%plot(1:nt,RMSE_noDA,'k-.',1:nt,RMSE_noDA_o,'r-.',1:nt,RMSE_noDA_n,'b-.','linewidth',3)
%legend('RMSE (all)','RMSE (obs)', 'RMSE (unobs)','RMSE (noDA all)', 'RMSE (noDA obs)','RMSE (noDA unobs)','fontsize',13,'location','northwest')
% legend('spread','spread (noDA)','RMSE','RMSE (noDA)','location','northwest')
%}


set(gca,'fontsize',24)
grid on
axis([0 nt 0 6])
% axis([200 500 1.5 7])
xlabel('timestep')
%title('PFF $$ H({\bf{x}}) = x_{(4a)} $$ $$\gamma = 1.25$$ ($$\alpha=0.01$$)','interpreter','latex','fontsize',30)
title('etkf realization 6 (model space)')
%title('LETKF $$ H({\bf{x}}) = x_{(4a)} $$ $$\gamma = 1.25$$','interpreter','latex','fontsize',30)
%title('G-MPF ($$\alpha = 60$$) $$ H({\bf{x}}) = x_{(4a)} $$','interpreter','latex','fontsize',35)
%title('MPF $$ H(x) = x^2 $$ obs freq = 20, G prior, variabe kernel','interpreter','latex','fontsize',22)
%title('MPF $$ H(x) = x_{4i-1}x_{4i} $$ obs freq = 20, G prior, fixed kernel','interpreter','latex','fontsize',22)
hold off