function peaks=FindPeak(data,epsilon)
global N_cluster 
%% == the first n columns of data contain the n-dimensional coordinates, the last column contains the value assigned to that point. 
[m,n]=size(data);
[~,ind]=sort(data(:,n));
peaks(1)=ind(1);
for i = 2:m
    iind=ind(i);
    dist=sqrt(sum((repmat(data(iind,1:n-1),length(peaks),1)-data(peaks,1:n-1)).^2,2));
    if (isempty(find(dist<epsilon)))
        peaks=[peaks,iind];
    else
        % disp('in the region')
    end
    if length(peaks)==N_cluster
        break;
    end
end



