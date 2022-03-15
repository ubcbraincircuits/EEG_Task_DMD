% This code plots the results

clear; clc; close all;

addpath('/Users/abhijitc/Documents/Software/eeglab2021.0');

data_sets = {'shamhceeg','shampd1eeg','shampd2eeg',...
            'stim7hceeg','stim7pd1eeg','stim7pd2eeg',...
            'stim8hceeg','stim8pd1eeg','stim8pd2eeg'};

nSub = [22,20,20,22,20,20,22,20,20];
FreqStr = {'Theta','Alpha','Beta','Gamma'};
fs = 1000;

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4];
horizon = Time_window';

DataType = 'rest';
RVal = 'ROpt';
NormVal = 'no_norm';
load(sprintf('SummaryData_%s_%s_%s.mat',DataType,RVal,NormVal));
plot_path = sprintf('Plots_2022_01_20/%s_data/Main_%s_%s/',DataType,RVal,NormVal);
mkdir(plot_path);
TestErrCh = TestErr;
TestErrCh(TestErrCh == 0) = NaN;

%% Plot Error - error norm

TestErr = squeeze(nanmean(nanmean(TestErrCh(:,1,:,:,:,:,:),1),4)); % error no norm
TDE_idx = [1,1,1,1];

Sub_ty_color = {'g','r','b'};
Stim_ty_str = {'sham','stim'};

%% PLot Test Error

fig = figure;
set(fig,'Position',[0 0 1440 900]);
error_threshold = NaN(1,length(FreqStr));

for stim_ty = 1:2
    for f = 1:length(FreqStr)
        s(f,stim_ty) = subplot(2,4,4*(stim_ty-1)+f);
        
        for sub_ty = 1:3
            ds_no = 2*(sub_ty-1)+stim_ty;
            N = sum(~isnan(squeeze(mean(TestErr(:,f,:,ds_no),1))));
            mean_plot = squeeze(nanmean(TestErr(:,f,:,ds_no),3)); 
            ste_plot = squeeze(nanstd(TestErr(:,f,:,ds_no)/sqrt(N),[],3));
            g = errorbar(horizon, mean_plot, ste_plot);
            set(g,'Color',Sub_ty_color{sub_ty},'LineStyle','none','LineWidth',2);
            hold on;
            fitObj = fit(horizon,mean_plot,'smoothingspline'); 
            h(sub_ty) = plot(horizon, fitObj(horizon));
            set(h(sub_ty),'Color',Sub_ty_color{sub_ty},'LineStyle','-','LineWidth',2);
        end
        
        temp = TestErr(:,f,:,:);
        error_threshold(f) = nanmean(temp,'all') + nanstd(temp(:));
        plot(horizon, ones(1,length(horizon)).*error_threshold(f),'k--','LineWidth',3);
        
        if stim_ty == 2
            xlabel('Time Horizon (No. of cycles)');
        elseif stim_ty == 1
            title(FreqStr{f});
        end
        if f == 1
            ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
        elseif f == 4 && stim_ty == 2
            legend(h,{'HC','PD, Med-off','PD, med-on'},'Location','Best');
        end
        
        %         ylim([0 0.3]);
    end
end

linkaxes(s(:));
sgtitle(sprintf('Test Error: data = %s, Error = div by (max-min), R = %s',strrep(NormVal,'_',' '),RVal));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/TestError_error_Ry.png',plot_path));

%% PLot Test Error 2

fig = figure;
set(fig,'Position',[0 0 1440 900]);
Sub_ty_str = {'HC','PD MedOff','PD MedOn'};
Stim_ty_color = {'k','m','c'};

for sub_ty = 1:3
    for f = 1:4
        s(f,sub_ty) = subplot(3,4,4*(sub_ty-1)+f);
        
        for stim_ty = 1:2
            ds_no = 2*(sub_ty-1)+stim_ty;
            N = sum(~isnan(squeeze(mean(TestErr(:,f,:,ds_no),1))));
            mean_plot = squeeze(nanmean(TestErr(:,f,:,ds_no),3)); 
            ste_plot = squeeze(nanstd(TestErr(:,f,:,ds_no)/sqrt(N),[],3));
            g = errorbar(horizon, mean_plot, ste_plot);
            set(g,'Color',Stim_ty_color{stim_ty},'LineStyle','none','LineWidth',2);
            hold on;
            fitObj = fit(horizon,mean_plot,'smoothingspline'); 
            h(stim_ty) = plot(horizon, fitObj(horizon));
            set(h(stim_ty),'Color',Stim_ty_color{stim_ty},'LineStyle','-','LineWidth',2);
        end
        
        plot(horizon, ones(1,length(horizon)).*error_threshold(f),'k--','LineWidth',3);
        
        if sub_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif sub_ty == 1
            title(FreqStr{f});
        end
        if f == 1
            ylabel(sprintf('%s',Sub_ty_str{sub_ty}));
        elseif f == 4 && sub_ty == 3
            legend(h,{'Sham','Stim'},'Location','Best');
        end
        
        %         ylim([0 0.3]);
    end
end

linkaxes(s(:));
sgtitle(sprintf('Test Error: data = %s, Error = div by (max-min), R = %s',strrep(NormVal,'_',' '),RVal));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/TestError_error_Ry_2.png',plot_path));

%% Get prediction horizon

TestErrMat = squeeze(TestErr); % Time x Freq x Sub x DS x TDE

for f = 1:4
    for ds = 1:6
        for sub = 1:22
            sub_data = TestErrMat(:,f,sub,ds);
            try
                fitObj = fit(horizon,sub_data,'smoothingspline');
                horizonFine = 0:0.001:4;
                yFine = fitObj(horizonFine);
                [~, idx] = min(abs(yFine - error_threshold(f)));
                PredHorizon(sub,ds,f) = horizonFine(idx);
            catch
                PredHorizon(sub,ds,f) = NaN;
            end
        end
    end
end

%% Plot correlation of pred hor. vs UPDRS

updrs = readtable('GVSinPD_UPDRS.xlsx');
UPDRS.total = updrs.UPDRSTotalScore;
UPDRS.rigidity = nansum([updrs.x3_3a,updrs.x3_3b,updrs.x3_3c,updrs.x3_3d,updrs.x3_3e],2);
UPDRS.tremor = nansum([updrs.x3_15a,updrs.x3_15b,updrs.x3_16a,updrs.x3_16b,updrs.x3_17a,updrs.x3_17b,updrs.x3_17c,updrs.x3_17d,updrs.x3_17e],2);
UPDRS.bradykinesia = updrs.x3_14;

PD_str = {'Med-Off','Med-On'};
updrs_type = {'total','rigidity','tremor','bradykinesia'};

for ty = 1:4
    for pd_ty = 1:2
        fig = figure;
        set(fig,'Position',[0 0 1440 900]);
        for f = 1:4
            for stim_ty = 1:2
                subplot(2,4,4*(stim_ty-1)+f);
                
                ds_no = 2*(pd_ty)+stim_ty;
                
                data = PredHorizon(:,ds_no,f);
                
                eval(sprintf('UPDRS_data = UPDRS.%s',updrs_type{ty}));
                
                scatter(data(1:20), UPDRS_data);
                [r,p] = corr(data(1:20), UPDRS_data);
                %axes square;
                if stim_ty == 2
                    xlabel('Time Horizon (No. of cycles)');
                end
                if f == 1
                    ylabel(sprintf('%s UPDRS',Stim_ty_str{stim_ty}));
                end
                title(sprintf('%s (r=%.2f,p=%.3f)',FreqStr{f},r,p));
                
            end
        end
        sgtitle(sprintf('UPDRS %s vs Prediction horizon for error threshold = %.2f, PD %s',updrs_type{ty},error_threshold(f),PD_str{pd_ty}));
        set(findall(gcf,'-property','FontSize'),'FontSize',15);
        saveas(fig,sprintf('%s/UPDRS_%s_Pred_horizon_thres_%.2f_%s.png',plot_path,updrs_type{ty},error_threshold(f),PD_str{pd_ty}));
    end
end

%% Plot bar with error

fig = figure;
set(fig,'Position',[0 0 1440 900]);
for f = 1:4
    for stim_ty = 1:2
        subplot(2,4,4*(stim_ty-1)+f);
        sub_ty = 1:3;
        ds_no = 2*(sub_ty-1)+stim_ty;
        
        data = PredHorizon(:,:,f);
        h = boxplot(data(:,ds_no));
        xticklabels({'HC','PD1','PD2'});
        
        if stim_ty == 2
            xlabel('Time Horizon (No. of cycles)');
        elseif stim_ty == 1
            title(sprintf('%s, thresh = %.2f',FreqStr{f},error_threshold(f)));
        end
        if f == 1
            ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
        end
        
    end
end
sgtitle(sprintf('Prediction horizon'));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/Pred_horizon.png',plot_path));

%% Plot bar with err

ColorStr = {'g','b','r'};
gap = 0.2222;

fig = figure('Position',[0 0 1440 900],'Color','white');
f = 3; stim_ty = 1;
sub_ty = 1:3;
ds_no = 2*(sub_ty-1)+stim_ty;
data = PredHorizon(:,:,f);
data_temp = data(:,ds_no);
n_temp = sum(~isnan(data_temp));

err_val = std(data_temp,'omitnan')./n_temp;
mean_val = mean(data_temp,'omitnan');

data = PredHorizon(:,:,f);
[b,e] = barwitherr([err_val; zeros(1,3)], [mean_val; zeros(1,3)], 'BarWidth', 0.5, 'FaceColor','w');
xlim([0.5 1.5]);
set(b,'LineWidth',5); set(e,'LineWidth',5);
for ii = 1:length(b)
    b(ii).EdgeColor = ColorStr{ii};
    e(ii).Color = ColorStr{ii};
%     if p_val(comp,ii) <= 0.05
%         
%     end
end
sigstar(x_groups{ii},p_val(comp,ii));
legend(b,{'HC','PD, MED-OFF','PD, MED-ON'});
box off;
h = gca;
ylim([2.4 3]);
h.XAxis.Visible = 'off';
title(sprintf('Prediction horizon (Beta, Sham)'));
set(findall(gcf,'-property','FontSize'),'FontSize',25);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
saveas(fig,sprintf('%s/Pred_horizon_beta_sham.png',plot_path));

%% Plot bar with error 2

fig = figure;
set(fig,'Position',[0 0 1440 900]);
for f = 1:4
    for sub_ty = 1:3
        subplot(3,4,4*(sub_ty-1)+f);
        stim_ty = 1:2;
        ds_no = 2*(sub_ty-1)+stim_ty;
        
        data = PredHorizon(:,:,f);
        h = boxplot(data(:,ds_no));
        xticklabels({'Sham','Stim'});
        
        if sub_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif sub_ty == 1
            title(sprintf('%s, thresh = %.2f',FreqStr{f},error_threshold(f)));
        end
        if f == 1
            ylabel(sprintf('%s',Sub_ty_str{sub_ty}));
        end
        
    end
end
suptitle(sprintf('Prediction horizon'));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/Pred_horizon2.png',plot_path));

%% ANOVA

ANOVA_path = sprintf('%s/ANOVA_plots/',plot_path);
mkdir(ANOVA_path);

for f = 1:4
    
    data = PredHorizon(:,:,f);
    sz = size(data);
    data = reshape(data,[sz(1:end-1),2,3]); % Sub x Patient Ty x Stim Ty
    [~, fig,~,f1,~] = cog_anova(data,[],{2,3,[2,3]},{'Sub #','Stim Type','Patient Type'});
    saveas(fig, sprintf('%s/ANOVA_%s_%.2f.png',ANOVA_path,FreqStr{f},error_threshold(f)));
    saveas(f1(1), sprintf('%s/Mult_comp_pat_%s_%.2f.png',ANOVA_path,FreqStr{f},error_threshold(f)));
    saveas(f1(2), sprintf('%s/Mult_comp_stim_%s_%.2f.png',ANOVA_path,FreqStr{f},error_threshold(f)));
    close all hidden;
    
end