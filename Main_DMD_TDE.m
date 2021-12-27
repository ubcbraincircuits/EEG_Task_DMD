%% This code runs DMD on all the subjects (usually run on GPU server)

clear; clc; close all;

fs = 500; % Sampling Freq
dt = 1/fs;

temp = who;
for ii = 1:length(temp)
    eval(sprintf('Params.%s = %s;', temp{ii},temp{ii}));
end

r_vals = [7,12,22,26]; %
TDE_vals = 1;

% Load data - make sure to have these in the same directory
load('READY4DL_500.mat');

data_sets = {'ALL_HC_GVS_OFF','ALL_HC_GVS_ON','ALL_PD_GVSOFF_MEDOFF',...
    'ALL_PD_GVSOFF_MEDON','ALL_PD_GVSON_MEDOFF','ALL_PD_GVSON_MEDON'};

mid_freq = [5.5, 11, 22, 39.5]; % Mid-frequencies for freq bands

% Design Filters for filtering. This takes some time, so I save the filters in .mat file 

try
    load('FilterParam.mat','FilterParam');
catch
    % Theta
    FilterParam(1) = designfilt('bandpassfir', 'StopbandFrequency1', 3, 'PassbandFrequency1', 4, 'PassbandFrequency2', 7, ...
        'StopbandFrequency2', 8, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', fs);
    % Alpha
    FilterParam(2) = designfilt('bandpassfir', 'StopbandFrequency1', 7, 'PassbandFrequency1', 8, 'PassbandFrequency2', 14, ...
        'StopbandFrequency2', 15, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', fs);
    % Beta
    FilterParam(3) = designfilt('bandpassfir', 'StopbandFrequency1', 14, 'PassbandFrequency1', 15, 'PassbandFrequency2', 29, ...
        'StopbandFrequency2', 30, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', fs);
    % Gamma
    FilterParam(4) = designfilt('bandpassfir', 'StopbandFrequency1', 29, 'PassbandFrequency1', 30, 'PassbandFrequency2', 49, ...
        'StopbandFrequency2', 50, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', fs);
    save('FilterParam.mat','FilterParam');
end

FreqStr = {'Theta','Alpha','Beta','Gamma'};

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4]; % Time windows to estimate errors

for ds = 1:length(data_sets)
    
    eval(sprintf('nSub = size(%s,1);',data_sets{ds}));
    eval(sprintf('dataset_val = %s;',data_sets{ds}));
    for ii = 1:size(dataset_val,2)/2
        for jj = 1:size(dataset_val,1)
            dataset_val2{jj,ii} = cat(1,dataset_val{jj,(ii-1)*2+1:2*ii});
        end
    end
    
    for sub = 1:nSub
        
        data = cat(3,dataset_val2{sub,:}); %s{sub};',data_sets{ds})); % data = 27 x 2000 x 10
        data = permute(data, [2,1,3]);
        data = cat(3,data(:,:,1:3),data,data(:,:,end-2:end)); % concatinate the 1st trial and last trial only to filter the data to avoid filter artifacts
        
        if any(isnan(data(:)))
            continue;
        end
        
        for f = 1:length(FreqStr)
            
            fprintf('%s Sub %d Freq %s started... \n',data_sets{ds},sub,FreqStr{f});
            [nCh, nSamples, nTrials] = size(data);
            dataMat = reshape(data, [nCh, nSamples*nTrials]);
            
            dataFilt = filtfilt(FilterParam(f), dataMat')'; % Filter data
            dataFilt = reshape(dataFilt, [nCh, nSamples, nTrials]);
            dataFilt(:,:,[1:3,end-2:end]) = []; % Reject 1-3 and last 3 trials due to filter artifcats
            
            [nCh, nSamples, nTrials] = size(dataFilt);
            
            save_path = sprintf('Results_DMD_test_rOpt_no_norm/%s/%s/Sub%d/',data_sets{ds},FreqStr{f},sub);
            mkdir(save_path);
            
            Params.nCh = nCh;
            Params.nSamples = nSamples;
            Params.nTrials = nTrials;
            Params.only_imag = true;
            Params.r = r_vals(f);
            Params.TDE = NaN;
            
            for tde = 1:length(TDE_vals)
                Params.TDE = round(TDE_vals(tde)*fs/mid_freq(f));
                parfor ii = 1:size(dataFilt,3)
                    Results_temp(ii) = mainFunction(dataFilt(:,:,ii), Params, 1, round(Time_window.*(fs./mid_freq(f))));
                end
                Results(tde,:) = Results_temp;
                clear Results_temp;
            end
            
            Results = Results';
            save([save_path,'Results.mat'],'Results');
            clear Results dataFilt;
            
        end
        clear data;
    end
end


%% Functions

function Results = mainFunction(data, Params, TrainTime, Time_window)

T = TrainTime*Params.fs; % No. of samples for training data

dataTrain = data(:,1:T);
dataTest = data(:,T-Params.TDE+1:end);
Ry2 = max(dataTest,[],2) - min(dataTest,[],2);

[S, f] = get_fft(dataTest',Params.fs,Params.fs);
[~, idx] = max(abs(S),[],1);
for ch = 1:size(S,2)
    Results.PhaseVal(ch) = angle(S(idx(ch),ch));
end

Results.Params = Params;

% Generate time delay embedding
Xaug_train = genTimeShiftEmbedding(dataTrain, Params.TDE);
Xaug_test = genTimeShiftEmbedding(dataTest, Params.TDE);

for ii = 1:length(Time_window)
    Xaug_train_cell{ii} = Xaug_train(:,1:Time_window(ii)+1);
    Xaug_test_cell{ii} = Xaug_test(:,1:Time_window(ii)+1);
end

Xaug_train_T1 = Xaug_train(:,1:end-1);
Xaug_train_T2 = Xaug_train(:,2:end);

Results.DmdStruct = runLowRankDMD(Xaug_train_T1, Xaug_train_T2, Params.r, Params.dt);

% Consider only imagine values
if Params.only_imag
    Results.DmdStruct.omega = 1i*imag(Results.DmdStruct.omega);
end

% Save the reconstruction errors
for ii = 1:length(Xaug_train_cell)
    [~, Results.reconErrorTrain(:,ii), ~, ~, ~] = predictDMD(Results.DmdStruct, Xaug_train_cell{ii}, Params.TDE, Ry2);
    [~, Results.reconErrorTest(:,ii), ~, ~, Results.reconErrorTestCh(:,:,ii)] = predictDMD(Results.DmdStruct, Xaug_test_cell{ii}, Params.TDE, Ry2);
end

end

function DmdStruct = runLowRankDMD(X,Y, r, dt)

[U, S, V] = svd(X, 'econ');
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);

Atilde = Ur'*Y*Vr/Sr;
[W, D] = eig(Atilde);
Phi = Y*Vr/Sr*W;

DmdStruct.S = S;
DmdStruct.r = r;
DmdStruct.dt = dt;
DmdStruct.A = Atilde;
DmdStruct.Phi = Phi;
DmdStruct.Lambda = diag(D);
DmdStruct.Freq = imag(log(diag(D))/dt)/(2*pi);
DmdStruct.omega = log(diag(D))/dt;

end

function Xaug = genTimeShiftEmbedding(X, nstacks)

% construct the augmented, shift-stacked data matrices
Xaug = [];
for st = 1:nstacks
    Xaug = [Xaug; X(:, st:end-nstacks+st)];
end

end

function [X, X2, X_big] = revTimeShiftEmbedding(Xaug, nstacks)

% reverse time delay embedding

ch = size(Xaug,1)/nstacks;
cols = size(Xaug, 2);
T = cols + nstacks - 1;

X_big = NaN(ch,T,cols);
X2 = [];

for col = 1:cols
    temp = reshape(Xaug(:,col),ch,nstacks);
    X_big(:,col:col+nstacks-1,col) = temp;
    if col == 1
        X2 = temp;
    else
        X2(:,end+1) = temp(:,end);
    end
end

X = nanmean(X_big, 3);

end

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
