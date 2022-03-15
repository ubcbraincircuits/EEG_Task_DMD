% This code runs through each folder and compliles necessary values

clear; clc; close all;
restoredefaultpath;
addpath("DMD_Functions/");

data_sets = {'shamhceeg','shampd1eeg','shampd2eeg',...
    'stim7hceeg','stim7pd1eeg','stim7pd2eeg',...
    'stim8hceeg','stim8pd1eeg','stim8pd2eeg'};

nSub = [22,20,20,22,20,20,22,20,20];
FreqStr = {'Theta','Alpha','Beta','Gamma'};
fs = 1000;

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4];
TDE_vals = 1;

RestTask = {'rest','task'};
RVals = {'rOpt'};
NormVals = {'noNorm'};

Freq_vals = [4 8;8 14;15 30;30 50];

for rt = 1:length(RestTask)
    for rv = 1:length(RVals)
        for nv = 1:length(NormVals)

            for ds = 1:length(data_sets)
                for f = 1:length(FreqStr)
                    for sub = 1:nSub(ds)
                        try
                            Res = load(sprintf('Results_DMD_%s_%s_%s/%s/%s/Sub%d/Results.mat',RestTask{rt},RVals{rv},NormVals{nv},data_sets{ds},FreqStr{f},sub));
                        catch
                            continue;
                        end
                        for tr = 1:length(Res.Results)
                            for tde = 1:length(TDE_vals)
                                TestErr(:,:,:,tr,f,sub,ds,tde) = Res.Results(tr,tde).reconErrorTestCh;
                                Omega_freq(:,tr,f,sub,ds,tde) = abs(Res.Results(tr,tde).DmdStruct.Freq); % modes x trial x freq band x sub x dataset x 1 TimeDelayEmbedding (TimeDelayEmbedding dimension is 1)
                                Phi(:,:,tr,f,sub,ds,tde) % EEG channels x modes x trial x freq band x sub x dataset x 1 (TimeDelayEmbedding dimension is 1)
                            end
                        end
                    end
                    disp(f);
                end
                clear Res;
                disp(ds);
            end

            save(sprintf('OmegaPhiData_%s_%s_%s.mat',RestTask{rt},RVals{rv},NormVals{nv}),"TestErr","Omega_freq","Phi"); %,'PhaseVals','Freq_Omega','Freq_Phi');
            clear TestErr Omega_freq Phi;
            fprintf('%s %s %s done \n\n',RestTask{rt},RVals{rv},NormVals{nv});

        end
    end
end

