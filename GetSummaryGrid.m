% This code runs through each folder and compliles necessary values

clear; clc; close all;

data_sets = {'shamhceeg','shampd1eeg','shampd2eeg',...
    'stim7hceeg','stim7pd1eeg','stim7pd2eeg',...
    'stim8hceeg','stim8pd1eeg','stim8pd2eeg'};

nSub = [22,20,20,22,20,20,22,20,20];
FreqStr = {'Theta','Alpha','Beta','Gamma'};
fs = 1000;
N = 5;
nDim = 50;

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4];
TDE_vals = 1;

nTr = 10;
TrainErr = NaN(length(Time_window),nTr,length(FreqStr),max(nSub),length(data_sets),length(TDE_vals));
TestErr = NaN(length(Time_window),nTr,length(FreqStr),max(nSub),length(data_sets),length(TDE_vals));
Freq = NaN(nDim,nTr,length(FreqStr),max(nSub),length(data_sets),length(TDE_vals));
EigVals = NaN(nDim,nTr,length(FreqStr),max(nSub),length(data_sets),length(TDE_vals));

for ds = 1:length(data_sets)
    for f = 1:length(FreqStr)
        for sub = 1:nSub(ds)
            try
                Res(sub,f) = load(sprintf('Results_DMD_test/%s/%s/Sub%d/Results.mat',data_sets{ds},FreqStr{f},sub));
            catch
                continue;
            end
            for tr = 1:nTr
                for tde = 1:length(TDE_vals)
                    TrainErr(:,tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).reconErrorTrain(5,:);
                    TestErr(:,tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).reconErrorTest(5,:);
                    Freq(1:length(Res(sub,f).Results(tr,tde).DmdStruct.Freq),tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).DmdStruct.Freq;
                    EigVals(1:length(Res(sub,f).Results(tr,tde).DmdStruct.Lambda),tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).DmdStruct.Lambda;
                end
            end
        end
        disp(f);
    end
    clear Res;
    disp(ds);
end

save(sprintf('SummaryData.mat'),'TrainErr','TestErr','Freq','EigVals');

