function peaks=FindPeak(data)
global coreNo epsilon
%% == the first n columns of data contain the n-dimensional coordinates, the last column contains the value assigned to that point. 
[m,n]=size(data);
[~,ind]=sort(data(:,n));
peaks(1)=ind(1);
index_counting=[ind(1)];
for i = 2:length(ind)
    iind=ind(i);
    dist=bsxfun(@minus,data(iind,1:3),data(peaks,1:3));
    dist=sqrt(sum(dist.^2,2));
    if (isempty(find(dist<epsilon)))
        peaks=[peaks,iind];
    end
    if length(peaks)==coreNo
        break;
    end
end



