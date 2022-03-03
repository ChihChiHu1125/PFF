% diagnostic for PFF
% RMSE and ensemble spread in observation space
% please run PFF_xxx.m first to generate the required information
% 2020/06/09
% 2022/03/01

%% Define input, parameters, and evaluate in obs space:

% define the state (trajectory/time series) of DA run and noDA run:
XX      = X;
XX_noDA = XnoDA;

% evaluate the state in observation space
title_name = 'RMSE obs space (linear sum)';

YY      = zeros(nt, ny_obs, np);  % DA ensemble in obs space
YY_noDA = zeros(nt, ny_obs, np);  % noDA ensemble in obs space 
YY_true = zeros(nt, ny_obs);      % truth in obs space

for i=1:nt
    for j=1:ny_obs
        inner_ind   = inner_domain{j};
        
        % modify the observation operator here:
        % linear_sum
        YY(i,j,:)      = H_linear_sum( XX(inner_ind,:,i), weights);
        YY_noDA(i,j,:) = H_linear_sum( XX_noDA(inner_ind,:,i), weights);
        YY_true(i,j)   = H_linear_sum( Xt(inner_ind, warm_nt+i), weights);

%         YY(i,j,:)      = H_square( XX(inner_ind,:,i));
%         YY_noDA(i,j,:) = H_square( XX_noDA(inner_ind,:,i));
%         YY_true(i,j)   = H_square( Xt(inner_ind, warm_nt+i));

    end
end

% RMSE
ens_mean      = squeeze(mean(YY,3));
ens_mean_noDA = squeeze(mean(YY_noDA,3));
RMSE_noDA     = zeros(nt,1);
RMSE          = zeros(nt,1);

for t=1:nt
    RMSE_noDA(t)   = sqrt(sum( (ens_mean_noDA(t,:)-YY_true(t,:)).^2 )/ny_obs);
    RMSE(t)        = sqrt(sum( (ens_mean(t,:)-YY_true(t,:)).^2 )/ny_obs);
end

% spread (mean spread across all obs)
spread      = zeros(nt,1);
spread_noDA = zeros(nt,1);

% calculate the variance first:
for t=1:nt
    spread(t)      = mean( std(YY(t,:,:),0,3) ); % standard deviation of the 3-rd dimension
    spread_noDA(t) = mean( std(YY_noDA(t,:,:),0,3) ); % standard deviation of the 3-rd dimension
end
%}


%% Plot section

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

plot(1:nt,spread,'color',[.8, .2, .2],'linewidth',2.5)
hold on
plot(1:nt,RMSE,'color',[.3 .3 .3],'linewidth',2.5)
plot(1:nt,spread_noDA,'-.','color',[.8, .2, .2],'linewidth',2.5)
plot(1:nt,RMSE_noDA,'-.','color',[.3, .3, .3],'linewidth',2.5)
legend('spread (DA)', 'RMSE (DA)','spread (noDA)', 'RMSE (noDA)','fontsize',13,'location','northwest')
set(gca,'fontsize',24)
grid on

% square
% axis([0 nt 0 100])
% axis([0 nt 0 70])

% exp
% axis([0 nt 0 10])

% abs or linear
axis([0 nt 0 3])

xlabel('timestep')
title(title_name,'interpreter','latex','fontsize',28)
hold off
%}