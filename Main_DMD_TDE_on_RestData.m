%% This code runs DMD on all the subjects (usually run on GPU server)

clear; clc; close all;
restoredefaultpath;
addpath("DMD_Functions/");

% Load data - make sure to have these in the same directory
load('/mnt/teamshare/MedData3/Wyatt/maryam_extract/hc_files.mat');
load('/mnt/teamshare/MedData3/Wyatt/maryam_extract/pdoff_files.mat');
load('/mnt/teamshare/MedData3/Wyatt/maryam_extract/pdon_files.mat');

shamhceeg = HCfiles(:,1,1);
shampd1eeg = PDOFFfiles(:,1,1);
shampd2eeg = PDONfiles(:,1,1);

stim7hceeg = HCfiles(:,2,1);
stim7pd1eeg = PDOFFfiles(:,2,1);
stim7pd2eeg = PDONfiles(:,2,1);

stim8hceeg = HCfiles(:,3,1);
stim8pd1eeg = PDOFFfiles(:,3,1);
stim8pd2eeg = PDONfiles(:,3,1);

for r_type = {'rOpt','r50'} % choose 'rOpt'or 'r50'
    for data_type = {'noNorm','zscore'} % choose 'noNorm' or 'zscore'

        clearvars -except r_type data_type *eeg;

        switch r_type{1}
            case 'rOpt'
                r_vals = [7,12,22,26]; %
            case 'r50'
                r_vals = [50,50,50,50]; %
        end

        Params.fs = 1000; % Sampling Freq
        Params.ds_factor = 1; % downsampling factor
        Params.dt = 1/Params.fs;
        Params.TDE_vals = 1;

        data_sets = {'shamhceeg','shampd1eeg','shampd2eeg',...
            'stim7hceeg','stim7pd1eeg','stim7pd2eeg',...
            'stim8hceeg','stim8pd1eeg','stim8pd2eeg'};

        Params.mid_freq = [5.5, 11, 22, 39.5]; % Mid-frequencies for freq bands

        % Design Filters for filtering. This takes some time, so I save the filters in .mat file

        try
            load('FilterParam.mat','FilterParam');
        catch
            % Theta
            FilterParam(1) = designfilt('bandpassfir', 'StopbandFrequency1', 3, 'PassbandFrequency1', 4, 'PassbandFrequency2', 7, ...
                'StopbandFrequency2', 8, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', Params.fs/Params.ds_factor);
            % Alpha
            FilterParam(2) = designfilt('bandpassfir', 'StopbandFrequency1', 7, 'PassbandFrequency1', 8, 'PassbandFrequency2', 14, ...
                'StopbandFrequency2', 15, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', Params.fs/Params.ds_factor);
            % Beta
            FilterParam(3) = designfilt('bandpassfir', 'StopbandFrequency1', 14, 'PassbandFrequency1', 15, 'PassbandFrequency2', 29, ...
                'StopbandFrequency2', 30, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', Params.fs/Params.ds_factor);
            % Gamma
            FilterParam(4) = designfilt('bandpassfir', 'StopbandFrequency1', 29, 'PassbandFrequency1', 30, 'PassbandFrequency2', 49, ...
                'StopbandFrequency2', 50, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', Params.fs/Params.ds_factor);
            save('FilterParam.mat','FilterParam');
        end

        FreqStr = {'Theta','Alpha','Beta','Gamma'};

        Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4]; % Time windows to estimate errors

        for ds = 1:length(data_sets)

            eval(sprintf('nSub = length(%s);',data_sets{ds}));

            for sub = 1:nSub

                eval(sprintf('dataStruct = %s{sub};',data_sets{ds})); 
                sz = size(dataStruct.data);
                
                data = reshape(dataStruct.data,[sz(1),2*Params.fs,sz(2)/(2*Params.fs)]); % data = 27 x 2000 x 30

                % Downsample data
                data = permute(data, [2,1,3]);
                data = downsample(data,Params.ds_factor);
                data = permute(data, [2,1,3]);
                Params.fs = Params.fs/Params.ds_factor;

                % save STD of data
                ResultsRawData.std = std(data,[],2);

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

                    switch data_type{1}
                        case 'zscore'
                            dataFilt = zscore(dataFilt,[],2);
                        case 'default'
                    end

                    % save STD of data
                    ResultsFiltData.std = std(dataFilt,[],2);

                    [nCh, nSamples, nTrials] = size(dataFilt);

                    save_path = sprintf('Results_DMD_rest_%s_%s/%s/%s/Sub%d/',r_type{1},data_type{1},data_sets{ds},FreqStr{f},sub);
                    mkdir(save_path);

                    Params.nCh = nCh;
                    Params.nSamples = nSamples;
                    Params.nTrials = nTrials;
                    Params.only_imag = true;
                    Params.r = r_vals(f);
                    Params.TDE = NaN;

                    for tde = 1:length(Params.TDE_vals)
                        Params.TDE = round(Params.TDE_vals(tde)*Params.fs/Params.mid_freq(f));
                        Time_window_f = round(Time_window.*(Params.fs./Params.mid_freq(f)));
                        parfor ii = 1:size(dataFilt,3)
                            Results_temp(ii) = mainDmdFunction(dataFilt(:,:,ii), Params, 1, Time_window_f);
                        end
                        Results(tde,:) = Results_temp;
                        clear Results_temp;
                    end

                    Results = Results';
                    save([save_path,'Results.mat'],'Results','ResultsFiltData','ResultsRawData','Params');
                    clear Results dataFilt;

                end
                clear data;
            end
        end

    end
end