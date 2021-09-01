% This code plots eigen values

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
plot_path = 'Plots_Eig/';
mkdir(plot_path);

%% Plot Eigen values

% EigVals (Dims x Tr x Freq x Sub x DS)

fig = figure;
set(fig,'Position',[0 0 1440 900]);

alp = 0.5;

for f = 1:4
    subplot(2,2,f);
    data = EigVals(:,:,f,:,1);
    s(1) = scatter(real(data(:)),imag(data(:)),'o','filled','MarkerFaceColor','g','MarkerFaceAlpha',alp);
    hold on;
    data = EigVals(:,:,f,:,2);
    s(2) = scatter(real(data(:)),imag(data(:)),'o','filled','MarkerFaceColor','r','MarkerFaceAlpha',alp);
    data = EigVals(:,:,f,:,3);
    s(3) = scatter(real(data(:)),imag(data(:)),'o','filled','MarkerFaceColor','b','MarkerFaceAlpha',alp);
    if f == 1
        legend(s,{'HC','PD med-off','PD med-on'});
    end
    xlabel('real'); ylabel('imag');
    title(FreqStr{f});
    %     xlim([-1.5 1.5]);
    %     ylim([-0.1,0.1]);
end
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(fig,sprintf('EigVals.png'));
