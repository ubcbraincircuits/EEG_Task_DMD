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
horizon = Time_window';

RVal = 'ROpt';
NormVal = 'no_norm';
load(sprintf('SummaryData_test_%s_%s.mat',RVal,NormVal));
plot_path = sprintf('Plots_2021_11_20/Main_%s_%s/Amp_mod/',RVal,NormVal);
mkdir(plot_path);

Freq_Omega(Freq_Omega == 0 | Freq_Phi == 0) = NaN;
Freq_Phi(Freq_Omega == 0 | Freq_Phi == 0) = NaN;

r_vals = [7,12,22,26];

Sub_ty_color = {'g','r','b'};
Stim_ty_str = {'sham','stim7','stim8'};
Sub_ty_str = {'HC','PD MedOff','PD MedOn'};
Stim_ty_color = {'k','m','c'};

%% plot results

fig = figure;
set(fig,'Position',[0 0 900 900]);
for f = 1:4
   subplot(2,2,f);
   f_omega = Freq_Omega(1:r_vals(f),:,:,f,:,:);
   f_phi = Freq_Phi(1:r_vals(f),:,:,f,:,:);
   scatter(f_omega(:), f_phi(:));
   axis equal;
   xlim([0 max([f_omega(:); f_phi(:)])]);
   ylim([0 max([f_omega(:); f_phi(:)])]);
   xlabel('Freq Omega');
   ylabel('Freq Phi');
   title(FreqStr{f});
end
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/freq.png',plot_path));

%%  Get freq assymmetry

for f = 1:4
    for ds = 1:9
        for sub = 1:22
            temp = abs(Freq_Omega(1:r_vals(f),1,:,f,sub,ds) - Freq_Phi(1:r_vals(f),1,:,f,sub,ds));
            freq_assym(sub,ds,f) = nanmean(temp(:));
        end
    end
end

freq_assym(freq_assym == 0) = NaN;

%% Plot bar with error

fig = figure;
set(fig,'Position',[0 0 1440 900]);
for f = 1:4
    for stim_ty = 1:3
        subplot(3,4,4*(stim_ty-1)+f);
        sub_ty = 1:3;
        ds_no = 3*(stim_ty-1)+sub_ty;
        
        data = freq_assym(:,:,f);
        h = boxplot(data(:,ds_no));
        xticklabels({'HC','PD1','PD2'});
        
        if stim_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif stim_ty == 1
            title(sprintf('%s',FreqStr{f}));
        end
        if f == 1
            ylabel(sprintf('%s',Stim_ty_str{stim_ty}));
        end
        
    end
end
suptitle(sprintf('Frequency asymmetry'));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/freq_asym.png',plot_path));

%% Plot bar with error 2

fig = figure;
set(fig,'Position',[0 0 1440 900]);
for f = 1:4
    for sub_ty = 1:3
        subplot(3,4,4*(sub_ty-1)+f);
        stim_ty = 1:3;
        ds_no = 3*(stim_ty-1)+sub_ty;
        
        data = freq_assym(:,:,f);
        h = boxplot(data(:,ds_no));
        xticklabels({'Sham','Stim-7','Stim-8'});
        
        if sub_ty == 3
            xlabel('Time Horizon (No. of cycles)');
        elseif sub_ty == 1
            title(sprintf('%s',FreqStr{f}));
        end
        if f == 1
            ylabel(sprintf('%s',Sub_ty_str{sub_ty}));
        end
        
    end
end
suptitle(sprintf('Frequency asymmetry'));
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('%s/freq_asym2.png',plot_path));

%% ANOVA

ANOVA_path = sprintf('%s/ANOVA_plots/',plot_path);
mkdir(ANOVA_path);

for f = 1:4
    
    data = freq_assym(:,:,f);
    sz = size(data);
    data = reshape(data,[sz(1:end-1),3,3]); % Sub x Patient Ty x Stim Ty
    [~, fig,~,f1,~] = cog_anova(data,[],{2,3,[2,3]},{'Sub #','Patient Type','Stim Type'});
    saveas(fig, sprintf('%s/ANOVA_%s.png',ANOVA_path,FreqStr{f}));
    saveas(f1(1), sprintf('%s/Mult_comp_pat_%s.png',ANOVA_path,FreqStr{f}));
    saveas(f1(2), sprintf('%s/Mult_comp_stim_%s.png',ANOVA_path,FreqStr{f}));
    close all hidden;
    
end

