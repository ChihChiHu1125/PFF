% Plot evolution of variable
XX = X;
var_id = 400;
da_intv = 20;
%name = ['no DA ensemble ($$x_{',num2str(var_id),'}$$)'];
name = ['noDA  $$H({\bf{x}}) = x^2_{(',num2str(var_id),')}$$'];
%name = ['LETKF $$H({\bf{x}}) = x^2_{(4a)}$$ ($$x_{',num2str(var_id),'}$$)'];
total_obs = nt/da_intv;
rng('default')

% for square observations
R        = 1*eye(dim/obs_den);
obs_rnd   = mvnrnd(zeros(ny_obs,1),R,total_obs)';
io_linear = 1; % if linear = 1, nonlinear = 0

% reproduce the observations:
observation = zeros(dim/obs_den, total_obs);
obs_time_ind = [da_intv:da_intv:nt];

if io_linear == 0
% square:
for t=1:total_obs
    %observation(:,t) = abs((Xt(obs_den:obs_den:end,warm_nt+da_intv*t))+obs_rnd(:,t));
    observation(:,t) = sqrt((Xt(obs_den:obs_den:end,warm_nt+da_intv*t).^2)+obs_rnd(:,t));
end

elseif io_linear == 1
% linear
for t=1:total_obs
    observation(:,t) = Xt(obs_den:obs_den:end,warm_nt+da_intv*t)+obs_rnd(:,t);
end
end


err = sqrt(0)*ones(total_obs,1); % the standard deviation of the observation error
figure;
set(gcf,'color','white')
set(gcf,'units','centimeters','position',[2 2 30 12])
%set(gcf,'units','centimeters','position',[2 2 40 12.5])


% for others:
plot(squeeze(XX(var_id,:,:))','linewidth',7,'color',[.5 .8 .5])
hold on
plot(1:nt, Xt(var_id, warm_nt+1:warm_nt+nt),'black','linewidth',4)
plot(squeeze(mean(XX(var_id,:,:),2)),'color',[.4 .4 .4],'linewidth',4)
if mod(var_id,4)==0
%{    
    for i=1:total_obs
        tmp  = observation(var_id/obs_den,i);
        tmp2 = Xt(var_id,warm_nt+da_intv*i);
        if tmp*tmp2>0
            errorbar(obs_time_ind(i), tmp, sqrt(0.5), ...
                's','markersize',8,'markeredgecolor','red','markerfacecolor','red')

            errorbar(obs_time_ind(i), -tmp, sqrt(0.5), ...
                's','markersize',8,'markeredgecolor','blue','markerfacecolor','blue')
        else
            errorbar(obs_time_ind(i), -tmp, sqrt(0.5), ...
                's','markersize',8,'markeredgecolor','red','markerfacecolor','red')

            errorbar(obs_time_ind(i), tmp, sqrt(0.5), ...
                's','markersize',8,'markeredgecolor','blue','markerfacecolor','blue')
        end
    end
%}
if io_linear ==1
    errorbar(obs_time_ind, observation(var_id/obs_den,:), err,'-s','markersize',8,'markeredgecolor','red','markerfacecolor','red')
elseif io_linear ==0
    for i=1:total_obs
        tmp  = observation(var_id/obs_den,i);
        tmp2 = Xt(var_id,warm_nt+da_intv*i);
        
        if tmp*tmp2>0
            plot(obs_time_ind(i), tmp, 's','color','red','linewidth',2)
            plot(obs_time_ind(i), -tmp,'s','color','blue','linewidth',2)
        else
            plot(obs_time_ind(i), -tmp, 's','color','red','linewidth',2)
            plot(obs_time_ind(i), tmp,'s','color','blue','linewidth',2)
        end
    end
    %{
    for i=1:total_obs
        tmp  = observation(var_id/obs_den,i);
        tmp2 = Xt(var_id,warm_nt+da_intv*i);
        if tmp*tmp2>0
            errorbar(obs_time_ind(i), tmp, sqrt(0.5)/abs(2*tmp2), ...
                '-s','markersize',8,'markeredgecolor','red','markerfacecolor','red','linewidth',1.5,'color','red')

            errorbar(obs_time_ind(i), -tmp, sqrt(0.5)/abs(2*tmp2), ...
                '-s','markersize',8,'markeredgecolor','blue','markerfacecolor','blue','linewidth',1.5,'color','blue')
        else
            errorbar(obs_time_ind(i), -tmp, sqrt(0.5)/abs(2*tmp2), ...
                '-s','markersize',8,'markeredgecolor','red','markerfacecolor','red','linewidth',1.5,'color','red')

            errorbar(obs_time_ind(i), tmp, sqrt(0.5)/abs(2*tmp2), ...
                '-s','markersize',8,'markeredgecolor','blue','markerfacecolor','blue','linewidth',1.5,'color','blue')
        end
    end
    %}
end

end

hold off
grid on
xticks([0:50:1500])
yticks([-20:4:20])
%axis([0 800 -22 18])
% axis([250 500 -12 12])
axis([0 300 -12 12])
% axis([0 40 -12 12])
xlabel('timestep')
set(gca,'fontsize',16)
%title('option3 (one of the observed variable)')
title(name,'interpreter','latex','fontsize',30)
%title('test11 $$ H(x) = x^2 $$ ($$x_4$$)','interpreter','latex','fontsize',16)