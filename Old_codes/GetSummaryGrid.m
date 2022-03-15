% This code runs through each folder and compliles necessary values

clear; clc; close all;

data_sets = {'ALL_HC_GVS_OFF','ALL_HC_GVS_ON','ALL_PD_GVSOFF_MEDOFF',...
    'ALL_PD_GVSOFF_MEDON','ALL_PD_GVSON_MEDOFF','ALL_PD_GVSON_MEDON'};

nSub = [22,20,20,22,20,20,22,20,20];
FreqStr = {'Theta','Alpha','Beta','Gamma'};
fs = 500;
N = 5;
nDim = 50;

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4];
TDE_vals = 1;

nTr = 10;
RVals = {'rOpt'};
NormVals = {'no_norm'};

for rv = 1:length(RVals)
    for nv = 1:length(NormVals)
        
        for ds = 1:length(data_sets)
            for f = 1:length(FreqStr)
                for sub = 1:nSub(ds)
                    try
                        Res(sub,f) = load(sprintf('Results_DMD_test_%s_%s/%s/%s/Sub%d/Results.mat',RVals{rv},NormVals{nv},data_sets{ds},FreqStr{f},sub));
                    catch
                        continue;
                    end
                    for tr = 1:nTr
                        for tde = 1:length(TDE_vals)
                            TestErr(:,:,:,tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).reconErrorTestCh;
                            PhaseVals(:,tr,f,sub,ds,tde) = Res(sub,f).Results(tr,tde).PhaseVal;
                            Phi = Res(sub,f).Results(tr,tde).DmdStruct.Phi;
                            for m = 1:size(Phi,2)
                                DMD_wavelet(:,:,m) = reshape(Phi(:,m),Res(sub,f).Results(tr,tde).Params.nCh,Res(sub,f).Results(tr,tde).Params.TDE);
                                %disp(m)
                                for ch = 1:Res(sub,f).Results(tr,tde).Params.nCh
                                    %disp(ch)
                                    dataR = real(DMD_wavelet(ch,:,m))';
                                    dataI = imag(DMD_wavelet(ch,:,m))';
                                    
                                    [S, freq] = get_fft(dataR.*hanning(length(dataR)), fs, 10*fs);
                                    [~, idx] = max(S);
                                    %disp(ch);
                                    Freq_Omega(m,ch,tr,f,sub,ds,tde) = abs(Res(sub,f).Results(tr,tde).DmdStruct.Freq(m));
                                    Freq_Phi(m,ch,tr,f,sub,ds,tde) = freq(idx);
                                    
                                end
                            end
                            clear DMD_wavelet;
                        end
                    end
                end
                disp(f);
            end
            clear Res;
            disp(ds);
        end
        
        save(sprintf('SummaryData_rest_%s_%s.mat',RVals{rv},NormVals{nv}),'TestErr','PhaseVals','Freq_Omega','Freq_Phi');
        clear TestErr PhaseVals Freq_Omega Freq_Phi;
        fprintf('%s %s done \n\n',RVals{rv},NormVals{nv});
        
    end
end

