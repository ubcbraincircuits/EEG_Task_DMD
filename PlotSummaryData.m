% This code plots the results

clear; clc; close all;

data_sets = {'shamhc','shampd1','shampd2',...
    'stim7hc','stim7pd1','stim7pd2',...
    'stim8hc','stim8pd1','stim8pd2'};
nSub = [22,20,20,22,20,20,22,20,20];
FreqStr = {'Theta','Alpha','Beta','Gamma'};
fs = 1000;
nDim = 50;
nTr = 10;

Time_window = [0.1:0.1:1 1.25:0.25:2 2.5:0.5:4];

load('SummaryData.mat');
plot_path = 'Plots/Main/';
mkdir(plot_path);

%% Plot Error

horizon = Time_window;

fig = figure;
set(fig,'Position',[0 0 1440 900]);

TrainErr = nanmean(TrainErr,2);
TestErr = nanmean(TestErr,2);
TDE_idx = [1,1,1,1];

Sub_ty_color = {'g','r','b'};
Stim_ty_str = {'sham','stim7','stim8'};

% Train Error

for stim_ty = 1:3
    for f = 1:4
        subplot(3,4,4*(stim_ty-1)+f);
        
        for sub_ty = 1:3
            ds_no = 3*(stim_ty-1)+sub_ty;
            N = sum(~isnan(squeeze(mean(TrainErr(:,:,f,:,ds_no,TDE_idx(f)),1))));
            h(sub_ty) = errorbar(horizon, squeeze(nanmean(TrainErr(:,:,f,:,ds_no,TDE_idx(f)),4)),...
                squeeze(nanstd(TrainErr(:,:,f,:,ds_no,TDE_idx(f))/sqrt(N),[],4)));
            set(h(sub_ty),'Color',Sub_ty_color{sub_ty},'LineStyle','-','LineWidth',2);
            hold on;
        end
        
        if stim_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif stim_ty == 1
            title(FreqStr{f});
        end
        if f == 1
            ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
        elseif f == 4 && stim_ty == 3
            legend(h,{'HC','PD1','PD2'},'Location','Best');
        end
        
%         ylim([0 0.3]);
    end
end

suptitle('Training Error');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/TrainError_plot_freq.png',plot_path));

% PLot Test Error

fig = figure;
set(fig,'Position',[0 0 1440 900]);

for stim_ty = 1:3
    for f = 1:4
        subplot(3,4,4*(stim_ty-1)+f);
        
        for sub_ty = 1:3
            ds_no = 3*(stim_ty-1)+sub_ty;
            N = sum(~isnan(squeeze(mean(TestErr(:,:,f,:,ds_no,TDE_idx(f)),1))));
            h(sub_ty) = errorbar(horizon, squeeze(nanmean(TestErr(:,:,f,:,ds_no,TDE_idx(f)),4)),...
                squeeze(nanstd(TestErr(:,:,f,:,ds_no,TDE_idx(f))/sqrt(N),[],4)));
            set(h(sub_ty),'Color',Sub_ty_color{sub_ty},'LineStyle','-','LineWidth',2);
            hold on;
        end
        
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

suptitle('Test Error');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/TestError_plot_freq.png',plot_path));

%% Get prediction horizon

TestErrMat = squeeze(TestErr); % Time x Freq x Sub x DS x TDE

ErrorThresh = [0.25,0.5,0.75, 1, 1.25, 1.5, 1.75, 2];

for f = 1:4
    for ds = 1:9
        for sub = 1:22
            sub_data = TestErrMat(:,f,sub,ds,TDE_idx(f));
            for th = 1:length(ErrorThresh)
                [~, idx] = min(abs(sub_data - ErrorThresh(th)));
                PredHorizon(sub,ds,f,th) = Time_window(idx);
            end
        end
    end
end

%% Plot bar with error 

for th = 1:length(ErrorThresh)
    fig = figure;
    set(fig,'Position',[0 0 1440 900]);
    for f = 1:4
        for stim_ty = 1:3
            subplot(3,4,4*(stim_ty-1)+f);
            sub_ty = 1:3; 
            ds_no = 3*(stim_ty-1)+sub_ty;
            
            data = PredHorizon(:,:,f,th);
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
    suptitle(sprintf('Prediction horizon for error threshold = %.2f',ErrorThresh(th)));
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    saveas(fig,sprintf('%s/Pred_horizon_thres_%.2f.png',plot_path,ErrorThresh(th)));
end

%% ANOVA

ANOVA_path = sprintf('%s/ANOVA_plots/',plot_path);
mkdir(ANOVA_path);

for th = 1:length(ErrorThresh)
    for f = 1:4
        
        data = PredHorizon(:,:,f,th);
        sz = size(data);
        data = reshape(data,[sz(1:end-1),3,3]); % Sub x Patient Ty x Stim Ty
        [~, fig,~,f1,~] = cog_anova(data,[],{2,3,[2,3]},{'Sub #','Patient Type','Stim Type'});
        saveas(fig, sprintf('%s/ANOVA_%s_%.2f.png',ANOVA_path,FreqStr{f},ErrorThresh(th)));
        saveas(f1(1), sprintf('%s/Mult_comp_pat_%s_%.2f.png',ANOVA_path,FreqStr{f},ErrorThresh(th)));
        saveas(f1(2), sprintf('%s/Mult_comp_stim_%s_%.2f.png',ANOVA_path,FreqStr{f},ErrorThresh(th)));
        close all hidden;
        
    end
end

%% Functions

function [p,fig,stats,f,p_mult] = cog_anova(data,group,effects,varnames)

anova_in = reshape(data,[],1);
size_val = size(data);

if ~exist('group','var') || isempty(group)
for i=1:ndims(data)
    group(:,i) = repmat(reshape(ones(prod(size_val(1:i-1)),1)*(1:size(data,i)),[],1),prod(size_val(i+1:end)),1);
end
end

group(isnan(anova_in),:) = [];
anova_in(isnan(anova_in),:) = [];

matrix = zeros(length(effects),size(group,2));

for i=1:length(effects)
matrix(i,effects{i}) = 1;
end

[p,~,stats,~,fig] = anovan_2(anova_in,group,'Random',1,'varnames',varnames,'model',matrix);

f(1) = figure;[p_mult] = multcompare(stats, 'Dimension',2);
f(2) = figure;[p_mult] = multcompare(stats, 'Dimension',3);

end
