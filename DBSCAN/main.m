%% Load Data
global epsm
[pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/160180408.bin');
% N_particle=dir([pathstr,'/', name, ext]).bytes/3/4; % 3: x,y,z-axis; 4:single precision
fid=fopen([pathstr,'/', name, ext]);
data=fread(fid,'single');
data=reshape(data,size(data,1)/3,3);
N_particle=size(data,1);
% figure,plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
subN=2000;%size(data,1);
ind=ceil(rand(subN,1)*N_particle);
X=data(ind,1:2);
figure, plot(data(:,1),data(:,2),'r.');
hold on; plot(X(:,1),X(:,2),'b*')

% X=[1 1; 1 2; 1 3; 1 4; 1 5; 2 3; 3 3; 4 3; 5 2; 5 3; 5 4; 6 3; 2 6; 3 5; 3 6; 3 7; 4 6];
epsm=0.003;
[ind_mbp,potential]=mbp(X,ones(subN,1),epsm);
%% Run DBSCAN Clustering Algorithm
clustering='dbscan';
if(strcmp(clustering,'kmeans'))
    IDX=kmeans(X,20);% 
    epsilon=0;
    MinPts=0;
elseif(strcmp(clustering,'dbscan'))
    epsilon=0.01;
    MinPts=10;
    IDX=dbscan(X,epsilon,MinPts);% 
end;
N_cluster=max(IDX);
fprintf('# of sub-clusters: %d \n',N_cluster);

h1=PlotClusterinResult(X, IDX,ind_mbp,epsilon,MinPts,clustering);
mbps=zeros(N_cluster,2);
center=zeros(N_cluster,2);
mass_sub=zeros(N_cluster,1);
for sub_cluster=1:N_cluster;
    ind_sub=find(IDX==sub_cluster);
    X_sub=X(ind_sub,:);
    subN=length(ind_sub);
    IDX_sub=dbscan(X_sub,epsilon,MinPts);
    ind_mbp_sub=mbp(X_sub,ones(subN,1),epsm);
    hold on; subplot(2,1,1);
    h2=plot(X_sub(ind_mbp_sub,1),X_sub(ind_mbp_sub,2),'bo');
    mbps(sub_cluster,:)=X_sub(ind_mbp_sub,:);
    mass_sub(sub_cluster)=subN;
    center(sub_cluster,:) = centroid(X_sub);
    % PlotClusterinResult(X, IDX,ind_mbp,epsilon,MinPts);
end
noise=X(IDX==0,:);
mbps=[mbps;noise];
mass_sub=[mass_sub;ones(size(noise,1),1)];
[ind_mbp_sub,potential_sub]=mbp(mbps,ones(size(mass_sub)),epsm);
 
points=unique([X(IDX==ind_mbp_sub,:);mbps],'rows');
[ind_2]=mbp(points,[ones(length(find(IDX==ind_mbp_sub)),1);mass_sub],epsm);

subplot(2,1,1);
hold on;
h2=plot(noise(:,1),noise(:,2),'bo');
h3=plot(mbps(ind_mbp_sub,1),mbps(ind_mbp_sub,2),'g*');
h4=plot(points(ind_2,1),points(ind_2,2),'m*');
hold off;
legend([h1,h2,h3,h4],'global-mbp','sub-mbp','one-mbp','Location', 'NorthEastOutside')





