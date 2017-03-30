%% Load Data
[pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/160180408.bin');
N_particle=str2num(name)/3/4; % 3: x,y,z-axis; 4:single precision
fid=fopen([pathstr,'/', name, ext]);
data=fread(fid,'single');
data=reshape(data,size(data,1)/3,3);
% figure,plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
subN=2000;
X=data(1:subN,1:2);

% X=[1 1; 1 2; 1 3; 1 4; 1 5; 2 3; 3 3; 4 3; 5 2; 5 3; 5 4; 6 3; 2 6; 3 5; 3 6; 3 7; 4 6];


%% Run DBSCAN Clustering Algorithm

epsilon=0.1;
MinPts=10;
% IDX=kmeans(X,10);% 
IDX=dbscan(X,epsilon,MinPts);% 
N_cluster=max(IDX);
fprintf('# of sub-clusters: %d \n',N_cluster);

[ind_mbp,potential]=mbp(X,ones(subN,1));

h1=PlotClusterinResult(X, IDX,ind_mbp,epsilon,MinPts);
% return;


mbps=zeros(N_cluster,2);
center=zeros(N_cluster,2);
mass_sub=zeros(N_cluster,1);
for sub_cluster=1:N_cluster;
    ind_sub=find(IDX==sub_cluster);
    X_sub=X(ind_sub,:);
    subN=length(ind_sub);
    IDX_sub=dbscan(X_sub,epsilon,MinPts);
    ind_mbp_sub=mbp(X_sub,ones(subN,1));
    hold on; subplot(2,1,1);
    h2=plot(X_sub(ind_mbp_sub,1),X_sub(ind_mbp_sub,2),'bo');
    mbps(sub_cluster,:)=X_sub(ind_mbp_sub,:);
    mass_sub(sub_cluster)=subN;
    center(sub_cluster,:) = centroid(X_sub);
    % PlotClusterinResult(X, IDX,ind_mbp,epsilon,MinPts);
end
[ind_mbp_sub,potential_sub]=mbp(mbps,mass_sub);
subplot(2,1,1);h3=plot(mbps(ind_mbp_sub,1),mbps(ind_mbp_sub,2),'g*');
legend([h1,h2,h3],'global-mbp','sub-mbp','one-mbp','Location', 'NorthEastOutside')





