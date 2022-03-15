function [Xdmd1, reconError, Xdmd2, X, reconError_Ch] = predictDMD(DmdStruct, X, nstacks, Ry2)

Phi = DmdStruct.Phi;
omega = DmdStruct.omega;
r = DmdStruct.r;
dt = DmdStruct.dt;

% Compute DMD mode amplitudes b
x1 = X(:,1);
b = Phi\x1;

% DMD reconstruction
mm1 = size(X, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1-1)*dt; % time vector

for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end

Xdmd = Phi * time_dynamics;

% Different kinds of reconstruction error

Ry = max(X,[],2) - min(X,[],2);
reconError(1) = mean(rms(X - real(Xdmd), 2)./Ry);

[~, X] = revTimeShiftEmbedding(X, nstacks);
[Xdmd1, Xdmd2] = revTimeShiftEmbedding(real(Xdmd), nstacks);

Ry = max(X,[],2) - min(X,[],2);
reconError(2) = mean(rms(X - Xdmd1, 2)./Ry);

% Normalized error using mean dmd signal
reconError_Ch(:,1) = rms(X(:,nstacks+1:end) - Xdmd1(:,nstacks+1:end),2)./Ry2;
reconError(3) = mean(reconError_Ch(:,1));

% Normalized error using last dmd signal
reconError_Ch(:,2) = rms(X(:,nstacks+1:end) - Xdmd2(:,nstacks+1:end),2)./Ry2;
reconError(4) = mean(reconError_Ch(:,2));

% This is the error using mean dmd signal (not normalized)
reconError_Ch(:,3) = rms(X(:,nstacks+1:end) - Xdmd1(:,nstacks+1:end),2);
reconError(5) = mean(reconError_Ch(:,3));

% This is the error using end dmd signal (not normalized)
reconError_Ch(:,4) = rms(X(:,nstacks+1:end) - Xdmd2(:,nstacks+1:end),2);
reconError(6) = mean(reconError_Ch(:,4));

end
