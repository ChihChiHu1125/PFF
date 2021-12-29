%% diagnostic (2020/04/24)
% rank histogram 
% 2020/06/21:
% rank historgram only for evaluating the forecast!
% prior: dim*np*total_obs
% rank histogram in the observation space!

% check the parameters are consistent:
warm_nt = 1000;
nt  = 1500;
dim = 1000;
np  =   20;

da_intv   = 20;              % how many timesteps between each observation to be assimilated
total_obs = nt/da_intv; % total number of observations (in time)

% observation variables:
all_ind = [1:1:dim];
obs_ind = [4:4:dim];
non_obs = setdiff(all_ind,obs_ind);
non_obs1 = [1:4:dim-3];
non_obs2 = [2:4:dim-2];
non_obs3 = [3:4:dim-1];

% test_ensemble name:
XX = abs_r1_prior;

%name1 = 'PFF obs var (obs space) $$ H({\bf{x}}) = |x_{(4a)}| $$';
%name1 = 'PFF obs var (obs space) $$ H({\bf{x}}) = e^{(\frac{x_{(4a)}}{6})} $$';
name1 = 'LETKF obs var (obs space) $$ H({\bf{x}}) = x_{(4a)}^2 $$';

div   = 30; % parameter for the plot (axis range)



% absolute value
%{
rk=zeros(dim,total_obs);
for t=1:total_obs
for i=1:dim
    tmp = sort([abs(squeeze(XX(i,:,t))),abs(Xt(i,warm_nt+da_intv*t))],'descend');
    rk(i,t) = find(tmp==abs(Xt(i,warm_nt+da_intv*t)));
end
end
%}

% exponential
%{
rk=zeros(dim,total_obs);
for t=1:total_obs
for i=1:dim
    tmp = sort([exp(squeeze(XX(i,:,t)/6)),exp(Xt(i,warm_nt+da_intv*t)/6)],'descend');
    rk(i,t) = find(tmp==exp(Xt(i,warm_nt+da_intv*t)/6));
end
end
%}

% square

rk=zeros(dim,total_obs);
for t=1:total_obs
for i=1:dim
    tmp = sort([squeeze(XX(i,:,t).^2),Xt(i,warm_nt+da_intv*t).^2],'descend');
    rk(i,t) = find(tmp==Xt(i,warm_nt+da_intv*t).^2);
end
end
%}



nbins = np + 1;
h = histogram(rk(obs_ind,:), nbins);
ct = h.Values;
%ct = ct/length(obs_ind);


set(gcf,'color','white') 
set(gcf,'units','centimeters','position',[5 5 20 14])
bar([0:1:20],ct/sum(ct)*100,'k')
%bar([0:1:20],etkf_sq_hist_obs/sum(etkf_sq_hist_obs)*100)
xticks([0:5:20])
xticklabels({'max','6','11','16','min'})
set(gca,'fontsize',22)
xlabel('rank','fontsize',26,'interpreter','latex')
ylabel('percentage','fontsize',26,'interpreter','latex')
title([name1],'interpreter','latex','fontsize',26)
axis([-1 21 0 div])




