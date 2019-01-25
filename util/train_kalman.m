function kf_model = train_kalman(X,Y,Z)
% x is input
% y is output
% z is measurement
if nargin < 3
    % measurement is previous state
    Z = Y(1:end-1,:);
    Y_new =  Y(2:end,:);
end

% silly geese doing it backwards
% X = X';
% Y = Y';
% Z = Z';
% X_new = X_new';

% system matrix
A = Y_new'*Z*inv(Z'*Z);
% observation matrix
C = X'*Y*(Y'*Y)^-1;
% model error covariance matrix
Q = (Y_new'*Y_new - A*(Z'*Y_new))/length(Y_new);
% observation covariance matrix
R = (X'*X - C*(Y'*X))/length(Y);



kf_model.A = A;
kf_model.C = C;
kf_model.Q = Q;
kf_model.R = R;





% % A
% [val1,val2] = deal(zeros(n_vars,n_vars));
% for t = 2:size(x,1)
%     val1 = val1 + x(t,:)'*x(t-1,:);
%     val2 = val2 + x(t-1,:)'*x(t-1,:);
% end
% A = val1 * inv(val2);
% 
% 
% % W
% [val1,val2] = deal(zeros(n_vars,n_vars));
% for t = 2:size(x,1)
%     val1 = val1 + x(t,:)'*x(t,:);
%     val2 = val2 + x(t-1,:)'*x(t,:);
% end
% W = (1/(N-1)) * (val1 - A * val2);
% 
% 
% % H
% [val1,val2] = deal(zeros(n_vars,n_vars));
% for t = 1:size(x,1)
%     val1 = val1 + z(t,:)'*x(t,:);
%     val2 = val2 + x(t,:)'*x(t,:);
% end
% H = val1 * inv(val2);
% 
% 
% % Q
% [val1,val2] = deal(zeros(n_vars,n_vars));
%     for t = 1:size(x,1)
%         val1 = val1 + z(t,:)'*z(t,:);
%         val2 = val2 + x(t,:)'*z(t,:);
%     end
% Q = (1/N) * (val1 - H * val2);
% 
% 
% 
% 
% 
% % Predict test trials
% for trial = 1:length(trial_data)
%     
%     % initial probability
%     Pk = eye(n_vars);
%     
%     x = trial_data(trial).(which_x);
%     z = trial_data(trial).(which_z);
%     y = trial_data(trial).(which_y);    
%     
%     xh = zeros(size(x));
%     xh(1,:) = y(1,:);
%     for t = 2:size(y,1)
%         Pk1 = Pk;
%         
%         % step 1
%         xhm = A * xh(t-1,:)';
%         Pkm = A * Pk1 * A' + W;
%         
%         % step 2
%         Kk = Pkm * H' * inv( H * Pkm * H' + Q);
%         Pk = (eye(size(Kk*H)) - Kk * H) * Pkm;
%         
%         xh(t,:) = xhm + Kk * (z(t,:)' - H * xhm);
%         
%     end
%     
%     trial_data(trial).(out_name) = xh;
%     
% end
% 
% 
