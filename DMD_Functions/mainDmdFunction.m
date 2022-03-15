function Results = mainDmdFunction(data, Params, TrainTime, Time_window)

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