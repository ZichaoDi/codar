%% Load Data
global epsm
rng('default')
[pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/160180408.bin');
% N_particle=dir([pathstr,'/', name, ext]).bytes/3/4; % 3: x,y,z-axis; 4:single precision
fid=fopen([pathstr,'/', name, ext]);
data=fread(fid,'single');
data=reshape(data,size(data,1)/3,3);
N_particle=size(data,1);
% figure,plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
subN=200;%size(data,1);
ind=ceil(rand(subN,1)*N_particle);
X=data(ind,1:2);
% figure, plot(data(:,1),data(:,2),'r.');
% hold on; plot(X(:,1),X(:,2),'b*')

% X=[1 1; 1 2; 1 3; 1 4; 1 5; 2 3; 3 3; 4 3; 5 2; 5 3; 5 4; 6 3; 2 6; 3 5; 3 6; 3 7; 4 6];
epsm=0.003;
tic;
[ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm);
toc;
%% Run DBSCAN Clustering Algorithm
tic;
clustering='kmeans';
if(strcmp(clustering,'kmeans'))
    N_cluster=3;
    IDX=kmeans(X,N_cluster);% 
    epsilon=0;
    MinPts=0;
elseif(strcmp(clustering,'dbscan'))
    epsilon=0.005;
    MinPts=10;
    IDX=dbscan(X,epsilon,MinPts);% 
    N_cluster=max(IDX)+1;
else
    epsilon=0.1;
    MinPts=0;
    N_cluster=10;
    ind_k=ceil(rand(N_cluster,1)*subN);
    IDX1=rangesearch(X,X(ind_k,:),epsilon);
    IDX=[];
    for j=1:N_cluster
        IDX(IDX1{j})=j;
    end
end;
h1=PlotClusterinResult(X,IDX,ind_mbp,epsilon,MinPts,clustering);
fprintf('# of sub-clusters: %d \n',N_cluster);
%%%====================== One step potential of particles in the same cluster
% mbps=zeros(N_cluster,2);
% center=zeros(N_cluster,2);
% super_potential=sparse(subN,1);
% mass_sp=sparse(N_cluster,1);
% for i=1:N_cluster
%     center(i,:) = centroid(X(IDX==i,:));
%     mass_sp(i)=length(find(IDX==i));
% end
% 
% for i=1:N_cluster;
%     ind_sub=find(IDX==i);
%     ind_other=setdiff([1:N_cluster],i);
%     X_sub=[X(ind_sub,:);center(ind_other,:)];
%     [~,potential]=mbp(X_sub,[ones(length(ind_sub),1);mass_sp(ind_other)],epsm);
%     super_potential(ind_sub)=potential(1:length(ind_sub));
% end
%%%============================================================
center=sparse(N_cluster,2);
mass_sp=ones(N_cluster,1);
super_potential=sparse(N_cluster,1);
for i=1:N_cluster
    center(i,:) = centroid(X(IDX==i,:));
    % mass_sp(i)=length(find(IDX==i));
end
for j=1:subN
    ind_other=[1:IDX(j)-1,IDX(j)+1:N_cluster];
    super_potential(j)=-sum(mass_sp(ind_other)./max(sqrt(sum(bsxfun(@minus,X(j,:),center(ind_other,:)).^2,2)),epsm));
end
%%%============================================================
[~,ind_super]=min(super_potential);
ind_super_old=ind_super;
%%%============================================================
super_potential_old=super_potential;
subplot(3,1,1);hold on; h2=plot(X(ind_super,1),X(ind_super,2),'go');
legend([h1,h2],'global-MBP','super-MBP');

[~,sortind]=sort(real_potential);
subplot(3,1,3), plot(real_potential(sortind),'r.-');hold on; plot(super_potential(sortind),'b.-');hold off;
legend('real-potential','super-potential')
%%%============================================================
n_level=1;
IDXtot=IDX;
ind_super_local=ind_super;
while(n_level<10)
    if(IDXtot(ind_super)==IDXtot(ind_super_old) & ind_mbp~=ind_super)
        disp(['Enter Phase ',num2str(n_level)]);
        Xsub=X(IDX==IDX(ind_super),:);
        subN=3;
        N_cluster=subN;
        mass_sp=ones(N_cluster,1);
        IDX_sub=find(IDX==IDX(ind_super));
        IDX=kmeans(Xsub,N_cluster);% 
        center=sparse(N_cluster,2);
        for i=1:N_cluster
            center(i,:) = centroid(Xsub(IDX==i,:));
            % mass_sp(i)=length(find(IDX==i));
        end
        super_potential=sparse(subN,1);
        for j=1:subN
            ind_other=[1:IDX(j)-1,IDX(j)+1:N_cluster];
            super_potential(IDX_sub(j))=-sum(mass_sp(ind_other)./max(sqrt(sum(bsxfun(@minus,Xsub(j,:),center(ind_other,:)).^2,2)),epsm));
        end
        [~,ind_super]=min(super_potential);
    elseif(ind_mbp==ind_super)
        disp('success')
    elseif(IDX(ind_super)~=IDX(ind_mbp))
        disp('diverge')
    end
    n_level=n_level+1;
end
subplot(3,1,1);hold on; plot(X(ind_super,1),X(ind_super,2),'bo');
figure,
plot(super_potential_old(sortind),'r.-');hold on; plot(super_potential(sortind),'b.-');hold off;
legend('real-potential','super-potential')
%%%======================== No particle, only super-particle
% for i=1:N_cluster
%     [~,potential]=mbp(X_sub,ones(subN,1),epsm);
%     super_potential(ind_sub)=potential;
%     hold on; subplot(2,1,1);
%     h2=plot(X_sub(ind_mbp_sub,1),X_sub(ind_mbp_sub,2),'bo');
%     mbps(i,:)=X_sub(ind_mbp_sub,:);
%     mass_sub(i)=subN;
%     % PlotClusterinResult(X, IDX,ind_mbp,epsilon,MinPts);
% end
% noise=X(IDX==0,:);
% mbps=[mbps;noise];
% mass_sub=[mass_sub;ones(size(noise,1),1)];
% [ind_mbp_sub,potential_sub]=mbp(mbps,ones(size(mass_sub)),epsm);
%  
% points=unique([X(IDX==ind_mbp_sub,:);mbps],'rows');
% [ind_2]=mbp(points,[ones(length(find(IDX==ind_mbp_sub)),1);mass_sub],epsm);
% 
% subplot(2,1,1);
% hold on;
% h2=plot(noise(:,1),noise(:,2),'bo');
% h3=plot(mbps(ind_mbp_sub,1),mbps(ind_mbp_sub,2),'g*');
% h4=plot(points(ind_2,1),points(ind_2,2),'m*');
% hold off;
% legend([h1,h2,h3,h4],'global-mbp','sub-mbp','one-mbp','Location', 'NorthEastOutside')
% 




