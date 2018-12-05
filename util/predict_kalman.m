function [STATEnew,Vnew,VVnew] = predict_kalman(kf_model,OBSERVATIONS,STATE0,V0)


A = kf_model.A;
C = kf_model.C;
Q = kf_model.Q;
R = kf_model.R;


STATEnew = zeros(size(OBSERVATIONS,1),size(STATE0,2));
[Vnew,VVnew] = deal(zeros(size(V0,1),size(V0,2),size(OBSERVATIONS,1)));

for t = 1:size(OBSERVATIONS,1)
    if t == 1 % define initial states for first time point
        STATEpred  = STATE0';
        Vpred = V0;
        Vprev = Vpred;
    else
        STATEprev = STATEnew(t-1,:)';
        Vprev = squeeze(Vnew(:,:,t-1));
        
        STATEpred = A * STATEprev;
        Vpred = A * Vprev * A' + Q;
    end
    
    err = OBSERVATIONS(t,:)' - C * STATEpred; % error (innovation)
    S = C * Vpred * C' + R;
    ss = size(V0,1);
    K = Vpred * C' * inv(S); % Kalman gain matrix
    
    % If there is no observation vector, set K = zeros(ss).
    STATEnew(t,:) = STATEpred + K * err;
    Vnew(:,:,t) = (eye(ss) - K * C) * Vpred;
    VVnew(:,:,t) = (eye(ss) - K * C) * A * Vprev;

    
    
end


