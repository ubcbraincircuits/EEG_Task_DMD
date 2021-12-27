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
f(2) = figure;[p_mult] = multcompare(stats, 'Dimension',[2,3]);

end