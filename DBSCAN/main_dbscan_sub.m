%% Load Data
global nrow epsm
rng('default')
% [pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/1069871097.bin');
[pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/771079093.bin');

% N_particle=dir([pathstr,'/', name, ext]).bytes/3/4; % 3: x,y,z-axis; 4:single precision
fid=fopen([pathstr,'/', name, ext]);
data=fread(fid,'single');
data=reshape(data,size(data,1)/3,3);
N_particle=size(data,1);
% figure,plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
subN=2000;%size(data,1);
totalN=subN;
ind=ceil(rand(subN,1)*N_particle);
X=data(ind,1:2);
epsm=0.003;
tic;
[ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm);
toc;
lb=min(real_potential);
ub=max(real_potential);
global_center=centroid(X);
figure,scatter(X(:,1),X(:,2),5,map1D(real_potential,[lb,ub]));
%% Run DBSCAN Clustering Algorithm
tic;
clustering='dbscan';
epsilon=0.002;
MinPts=10;
[IDX,c]=dbscan(X,epsilon,MinPts);% 
toc;
tic;
N_cluster=max(IDX);
mass_sp=sparse(N_cluster,1);
for i=1:N_cluster
    mass_sp(i)=length(find(IDX==i));
end
[~,ind_target]=max(mass_sp);
ind_sub=find(IDX==ind_target);
N_sp=mass_sp(ind_target);
super_potential=sparse(N_sp,1);
disp(['Cardinality of maximum super particle ',num2str(N_sp)]);
for i=1:N_sp
    ind_other=[1:i-1,i+1:subN];
    super_potential(i)=-sum(ones(length(ind_other),1)./max(sqrt(sum(bsxfun(@minus,X(ind_sub(i),:),X(ind_other,:)).^2,2)),epsm));
end
toc;
[~,ind_super]=min(super_potential);
nrow=3;
h1=PlotClusterinResult(X,IDX,ind_mbp,epsilon,MinPts,clustering);
fprintf('# of sub-clusters: %d \n',N_cluster);
figure(findobj('Tag','ml'));
subplot(nrow,1,1);hold on; h2=plot(X(ind_sub(ind_super),1),X(ind_sub(ind_super),2),'go');

% %%%============================================================
legend([h1,h2],'global-MBP','super-MBP');




