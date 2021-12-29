function X_out = L9640_RK4(X_in, dt,F)
% integrate for L96 model using RK4
% X_in is the input, dt is the time resolution
% X_out is the output

[dim, np] = size(X_in);
k1 = zeros(dim,np);
k2 = zeros(dim,np);
k3 = zeros(dim,np);
k4 = zeros(dim,np);
X_out = zeros(dim,np);

% parameter for the model
% F = 8;

tmp_b = X_in; % before integration
X_p1  = [tmp_b(2:end,:); tmp_b(1,:)];
X_00  = [tmp_b];
X_n1  = [tmp_b(end,:); tmp_b(1:end-1,:)];
X_n2  = [tmp_b(end-1:end,:); tmp_b(1:end-2,:)];
    
k1 = (X_p1-X_n2).*X_n1 - X_00 + F;
    
tmp_b = X_in + 0.5*k1*dt;
X_p1  = [tmp_b(2:end,:); tmp_b(1,:)];
X_00  = [tmp_b];
X_n1  = [tmp_b(end,:); tmp_b(1:end-1,:)];
X_n2  = [tmp_b(end-1:end,:); tmp_b(1:end-2,:)];
    
k2 = (X_p1-X_n2).*X_n1 - X_00 + F;
    
tmp_b = X_in + 0.5*k2*dt;
X_p1  = [tmp_b(2:end,:); tmp_b(1,:)];
X_00  = [tmp_b];
X_n1  = [tmp_b(end,:); tmp_b(1:end-1,:)];
X_n2  = [tmp_b(end-1:end,:); tmp_b(1:end-2,:)];
    
k3 = (X_p1-X_n2).*X_n1 - X_00 + F;    
    
tmp_b = X_in + k3*dt;
X_p1  = [tmp_b(2:end,:); tmp_b(1,:)];
X_00  = [tmp_b];
X_n1  = [tmp_b(end,:); tmp_b(1:end-1,:)];
X_n2  = [tmp_b(end-1:end,:); tmp_b(1:end-2,:)];
    
k4 = (X_p1-X_n2).*X_n1 - X_00 + F; 
    
X_out = X_in + (dt/6)*(k1+2*k2+2*k3+k4);
  

end