function [X, X2, X_big] = revTimeShiftEmbedding(Xaug, nstacks)

% reverse time delay embedding

ch = size(Xaug,1)/nstacks;
cols = size(Xaug, 2);
T = cols + nstacks - 1;

X_big = NaN(ch,T,cols);
X2 = [];

for col = 1:cols
    temp = reshape(Xaug(:,col),ch,nstacks);
    X_big(:,col:col+nstacks-1,col) = temp;
    if col == 1
        X2 = temp;
    else
        X2(:,end+1) = temp(:,end);
    end
end

X = nanmean(X_big, 3);
end