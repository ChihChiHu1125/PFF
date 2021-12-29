%% diagnostic (2020/04/24)
% rank histogram 
% 2020/05/27 revise:
% rank historgram only for evaluating the forecast!
% prior: dim*np*total_obs

% check the parameters are consistent:
warm_nt = 1000;
nt  = 200;
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
XX = X;
name1 = 'PFF $$ H({\bf{x}}) = x_{(4a)}$$ ($$\alpha = 1$$) (obs)';
name2 = 'PFF $$ H({\bf{x}}) = x_{(4a)}$$ ($$\alpha = 1$$) (non-obs';
%name1 = 'LETKF $$ H({\bf{x}}) = x_{(4a)}$$ $$\gamma = 1$$ (obs)';
%name2 = 'LETKF $$ H({\bf{x}}) = x_{(4a)}$$ $$\gamma = 1$$ (non-obs';
%name1 = 'MPF (B) ($$\alpha = 0.05$$) $$ H({\bf{x}}) = x_{(4a)}$$ (obs)';
%name2 = 'MPF (B) ($$\alpha = 0.05$$) $$ H({\bf{x}}) = x_{(4a)}$$ (non-obs';
%name1 = 'MPF $$ H({\bf{x}}) = x^2_{(4a)}$$ GM prior (obs freq = 5) ni=250 (obs)';
%name2 = 'MPF $$ H({\bf{x}}) = x^2_{(4a)}$$ GM prior (obs freq = 5) ni=250 (non-obs';
div   = 40; % parameter for the plot (axis range)

% Posterior:
%{
rk=zeros(dim,total_obs);
for t=1:total_obs
for i=1:dim
    tmp = sort([squeeze(XX(i,:,da_intv*t)),Xt(i,warm_nt+da_intv*t)],'descend');
    rk(i,t) = find(tmp==Xt(i,warm_nt+da_intv*t));
end
end
%}

% Prior:


rk=zeros(dim,total_obs);
for t=1:total_obs
for i=1:dim
    tmp = sort([squeeze(XX(i,:,t)),Xt(i,warm_nt+da_intv*t)],'descend');
    rk(i,t) = find(tmp==Xt(i,warm_nt+da_intv*t));
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


%{
h = histogram(rk(non_obs,:), nbins);
ct = h.Values;
ct = ct/length(non_obs);


set(gcf,'color','white') 
set(gcf,'units','centimeters','position',[5 5 20 14])
bar([0:1:20],ct)
xticks([0:5:20])
xticklabels({'max','6','11','16','min'})
set(gca,'fontsize',15)
xlabel('rank','fontsize',20,'interpreter','latex')
ylabel('count','fontsize',20,'interpreter','latex')
title('non-observed variables')
%}

%{
h = histogram(rk(non_obs1,:), nbins);
ct = h.Values;
%ct = ct/length(non_obs1);


set(gcf,'color','white') 
set(gcf,'units','centimeters','position',[5 5 20 14])
bar([0:1:20],ct/sum(ct)*100)
%bar([0:1:20],etkf_sq_hist_no_obs1/sum(etkf_sq_hist_no_obs1)*100)
xticks([0:5:20])
xticklabels({'max','6','11','16','min'})
set(gca,'fontsize',15)
xlabel('rank','fontsize',20,'interpreter','latex')
ylabel('percentage','fontsize',20,'interpreter','latex')
title([name2,' 1)'],'interpreter','latex','fontsize',25)
axis([-1 21 0 div])
pause


h = histogram(rk(non_obs2,:), nbins);
ct = h.Values;
%ct = ct/length(non_obs2);


set(gcf,'color','white') 
set(gcf,'units','centimeters','position',[5 5 20 14])
bar([0:1:20],ct/sum(ct)*100)
%bar([0:1:20],etkf_sq_hist_no_obs2/sum(etkf_sq_hist_no_obs2)*100)
xticks([0:5:20])
xticklabels({'max','6','11','16','min'})
set(gca,'fontsize',15)
xlabel('rank','fontsize',20,'interpreter','latex')
ylabel('percentage','fontsize',20,'interpreter','latex')
title([name2,' 2)'],'interpreter','latex','fontsize',25)
axis([-1 21 0 div])
pause


h = histogram(rk(non_obs3,:), nbins);
ct = h.Values;
%ct = ct/length(non_obs3);


set(gcf,'color','white') 
set(gcf,'units','centimeters','position',[5 5 20 14])
bar([0:1:20],ct/sum(ct)*100)
%bar([0:1:20],etkf_sq_hist_no_obs3/sum(etkf_sq_hist_no_obs3)*100)
xticks([0:5:20])
xticklabels({'max','6','11','16','min'})
set(gca,'fontsize',15)
xlabel('rank','fontsize',20,'interpreter','latex')
ylabel('percentage','fontsize',20,'interpreter','latex')
title([name2,' 3)'],'interpreter','latex','fontsize',25)
axis([-1 21 0 div])
%}
