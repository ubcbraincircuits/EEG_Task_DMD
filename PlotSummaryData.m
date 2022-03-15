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
NormVal = 'noNorm';
load(sprintf('SummaryData_%s_%s_%s.mat',DataType,RVal,NormVal));
plot_path = sprintf('Plots_2022_01_20/%s_data/Main_%s_%s/',DataType,RVal,NormVal);
mkdir(plot_path);
TestErrCh = TestErr;
TestErrCh(TestErrCh == 0) = NaN;

%% Plot Error - error norm

TestErr = squeeze(nanmean(nanmean(TestErrCh(:,1,:,1:10,:,:,:),1),4)); % error no norm
TDE_idx = [1,1,1,1];

Sub_ty_color = {'g','r','b'};
Stim_ty_str = {'sham','stim7','stim8'};

%% PLot Test Error

fig = figure;
set(fig,'Position',[0 0 1440 900]);

for stim_ty = 1:3
    for f = 1:4
        s(f,stim_ty) = subplot(3,4,4*(stim_ty-1)+f);

        for sub_ty = 1:3
            ds_no = 3*(stim_ty-1)+sub_ty;
            N = sum(~isnan(squeeze(mean(TestErr(:,f,:,ds_no),1))));
            h(sub_ty) = errorbar(horizon, squeeze(nanmean(TestErr(:,f,:,ds_no),3)),...
                squeeze(nanstd(TestErr(:,f,:,ds_no)/sqrt(N),[],3)));
            set(h(sub_ty),'Color',Sub_ty_color{sub_ty},'LineStyle','-','LineWidth',2);
            hold on;
        end

        temp = TestErr(:,f,:,:);
        error_threshold(f) = nanmean(temp,'all') + nanstd(temp(:));
        plot(horizon, ones(1,length(horizon)).*error_threshold(f),'k--','LineWidth',3);

        if stim_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif stim_ty == 1
            title(FreqStr{f});
        end
        if f == 1
            ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
        elseif f == 4 && stim_ty == 3
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

        for stim_ty = 1:3
            ds_no = 3*(stim_ty-1)+sub_ty;
            N = sum(~isnan(squeeze(mean(TestErr(:,f,:,ds_no),1))));
            h(stim_ty) = errorbar(horizon, squeeze(nanmean(TestErr(:,f,:,ds_no),3)),...
                squeeze(nanstd(TestErr(:,f,:,ds_no)/sqrt(N),[],3)));
            set(h(stim_ty),'Color',Stim_ty_color{stim_ty},'LineStyle','-','LineWidth',2);
            hold on;
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
            legend(h,{'Sham','Stim-7','Stim-8'},'Location','Best');
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
    for ds = 1:9
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

UPDRS = [54,42,15,47,42,40,66,62,21,35,31,41,40,54,22,37,38,30,29,29];
PD_str = {'Med-Off','Med-On'};

for pd_ty = 1:2
    fig = figure;
    set(fig,'Position',[0 0 1440 900]);
    for f = 1:4
        for stim_ty = 1:3
            subplot(3,4,4*(stim_ty-1)+f);

            ds_no = 3*(stim_ty-1)+pd_ty+1;

            data = PredHorizon(:,ds_no,f);

            scatter(data(1:20), UPDRS);
            [r,p] = corr(data(1:20), UPDRS');
            %axes square;
            if stim_ty == 3
                xlabel('Time Horizon (No. of cycles)');
            end
            if f == 1
                ylabel(sprintf('%s UPDRS',Stim_ty_str{stim_ty}));
            end
            title(sprintf('%s (r=%.2f,p=%.3f)',FreqStr{f},r,p));

        end
    end
    sgtitle(sprintf('UPDRS vs Prediction horizon for PD %s',PD_str{pd_ty}));
    set(findall(gcf,'-property','FontSize'),'FontSize',15);
    saveas(fig,sprintf('%s/UPDRS_Pred_horizon_%s.png',plot_path,PD_str{pd_ty}));
end

%% Plot bar with error

    fig = figure;
    set(fig,'Position',[0 0 1440 900]);
    for f = 1:4
        for stim_ty = 1:3
            subplot(3,4,4*(stim_ty-1)+f);
            sub_ty = 1:3;
            ds_no = 3*(stim_ty-1)+sub_ty;

            data = PredHorizon(:,:,f);
            h = boxplot(data(:,ds_no));
            xticklabels({'HC','PD1','PD2'});

            if stim_ty == 3
                xlabel('Time Horizon (No. of cycles)');
            elseif stim_ty == 1
                title(FreqStr{f});
            end
            if f == 1
                ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
            end

        end
    end
    sgtitle(sprintf('Prediction horizon'));
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(fig,sprintf('%s/Pred_horizon.png',plot_path));

%% Plot bar with error 2

fig = figure;
set(fig,'Position',[0 0 1440 900]);
for f = 1:4
    for sub_ty = 1:3
        subplot(3,4,4*(sub_ty-1)+f);
        stim_ty = 1:3;
        ds_no = 3*(stim_ty-1)+sub_ty;

        data = PredHorizon(:,:,f);
        h = boxplot(data(:,ds_no));
        xticklabels({'Sham','Stim-7','Stim-8'});

        if sub_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif sub_ty == 1
            title(FreqStr{f});
        end
        if f == 1
            ylabel(sprintf('%s',Sub_ty_str{sub_ty}));
        end

    end
end
sgtitle(sprintf('Prediction horizon'));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/Pred_horizon2.png',plot_path));

%% ANOVA

ANOVA_path = sprintf('%s/ANOVA_plots/',plot_path);
mkdir(ANOVA_path);

for f = 1:4

    data = PredHorizon(:,:,f);
    sz = size(data);
    data = reshape(data,[sz(1:end-1),3,3]); % Sub x Patient Ty x Stim Ty
    [~, fig,~,f1,~] = cog_anova(data,[],{2,3,[2,3]},{'Sub #','Patient Type','Stim Type'});
    saveas(fig, sprintf('%s/ANOVA_%s.png',ANOVA_path,FreqStr{f}));
    saveas(f1(1), sprintf('%s/Mult_comp_pat_%s.png',ANOVA_path,FreqStr{f}));
    saveas(f1(2), sprintf('%s/Mult_comp_stim_%s.png',ANOVA_path,FreqStr{f}));
    saveas(f1(3), sprintf('%s/Mult_comp_stim_x_pat_%s.png',ANOVA_path,FreqStr{f}));
    close all hidden;

end