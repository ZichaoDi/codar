%% Load Data
global nrow epsm
rng('default')
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
% figure, plot(data(:,1),data(:,2),'r.');
% hold on; plot(X(:,1),X(:,2),'b*')

% X=[1 1; 1 2; 1 3; 1 4; 1 5; 2 3; 3 3; 4 3; 5 2; 5 3; 5 4; 6 3; 2 6; 3 5; 3 6; 3 7; 4 6];
epsm=0.003;
[ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm);
lb=min(real_potential);
ub=max(real_potential);
global_center=centroid(X);
figure,scatter(X(:,1),X(:,2),5,map1D(real_potential,[lb,ub]));
%% Run DBSCAN Clustering Algorithm
clustering='dbscan';
if(strcmp(clustering,'dbscan'))
    epsilon=0.005;
    MinPts=10;
    [IDX,c]=dbscan(X,epsilon,MinPts);% 
    N_cluster=max(IDX);
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
nrow=4;
h1=PlotClusterinResult(X,IDX,ind_mbp,epsilon,MinPts,clustering);
fprintf('# of sub-clusters: %d \n',N_cluster);
mass_sp=ones(N_cluster,1);
super_potential=sparse(N_cluster,1);
    center=c;
    for i=1:N_cluster
        mass_sp(i)=1/length(find(IDX==i));
    end
%%==================================================
for j=1:subN
    if(N_cluster>1)
        ind_other=[1:IDX(j)-1,IDX(j)+1:N_cluster];
        super_potential(j)=-sum(mass_sp(ind_other)./max(sqrt(sum(bsxfun(@minus,X(j,:),center(ind_other,:)).^2,2)),epsm));
    else
        return;
    end
end
%%%============================================================
[~,ind_super]=min(super_potential);
ind_super_old=ind_super;
%%==================================================
f_dy=findobj('Tag','dynamics');
if(isempty(f_dy))
    f_dy=figure('Tag','dynamics');
else
    figure(f_dy)
end
max_level=2;
subplot(max_level,2,1);
Colors=hsv(N_cluster);
for jj=1:N_cluster
plot(X(IDX==jj,1),X(IDX==jj,2),'.','color',Colors(jj,:));
hold on;
end
plot(X(ind_mbp,1),X(ind_mbp,2),'k*',center(:,1),center(:,2),'go',X(ind_super,1),X(ind_super,2),'kd','markersize',7);% text(center(:,1),center(:,2),num2str(mass_sp));
subplot(max_level,2,2)
scatter(X(:,1),X(:,2),5,map1D(super_potential,[lb,ub]));
%%%============================================================
super_potential_old=super_potential;
figure(findobj('Tag','ml'));
subplot(nrow,1,1);hold on; h2=plot(X(ind_super,1),X(ind_super,2),'go');

[~,sortind]=sort(real_potential);
[~,sortind1]=sort(super_potential);
subplot(nrow,1,3), plot(real_potential(sortind),'r.-');hold on; plot(map1D(super_potential(sortind),[lb,ub]),'b.-');hold off;
legend('real-potential','super-potential')
subplot(nrow,1,4), plot(real_potential(sortind1),'r.-');hold on; plot(map1D(super_potential(sortind1),[lb,ub]),'b.-');hold off;
legend('real-potential','super-potential')
% %%%============================================================
n_level=1;
ind_track=[];
n_cluster=N_cluster;
while(n_level<max_level)
    if(ind_super_old~=ind_mbp)
        disp(['Enter Phase ',num2str(n_level)]);
        ind_sub=find(IDX==IDX(ind_super));
        Xsub=X(ind_sub,:);
        epsilon=epsilon;
        [IDX_sub,c]=dbscan(Xsub,epsilon,MinPts/2);% 
        N_cluster=max(IDX_sub);
        if(N_cluster==0 | IDX(ind_super)==0)
            disp('no further clustering')
            break;
        end
        center(IDX(ind_super),:)=c(1,:);
        IDX_old=IDX;
        IDX_temp=IDX;
        for jj=1:N_cluster
            if jj==1
                IDX_temp(ind_sub(IDX_sub==jj))=IDX(ind_super);
            else
                IDX_temp(ind_sub(IDX_sub==jj))=jj+n_cluster-1;
            end
        end
        IDX=IDX_temp;
            
        center=[center;c(2:end,:)];
        n_cluster=n_cluster-1+N_cluster;
        mass_sp=ones(n_cluster,1);
        for i=1:n_cluster
            mass_sp(i)=1/length(find(IDX==i));
        end
        for j=1:subN
            ind_other=[1:IDX(j)-1,IDX(j)+1:N_cluster];
            super_potential(j)=-sum(mass_sp(ind_other)./max(sqrt(sum(bsxfun(@minus,X(j,:),center(ind_other,:)).^2,2)),epsm));
        end
        [~,ind_super]=min(super_potential);
        if(ind_super_old==ind_super)
            disp('trapped')
            break;
        end
        ind_super_old=ind_super;
        %%==========================================
        figure(findobj('Tag','dynamics'));
        subplot(max_level,2,2*n_level+1),
        Colors=hsv(n_cluster);
        for jj=1:n_cluster
            plot(X(IDX==jj,1),X(IDX==jj,2),'.','color',Colors(jj,:));
            hold on;
        end
        plot(X(ind_mbp,1),X(ind_mbp,2),'k*',center(:,1),center(:,2),'go',X(ind_super,1),X(ind_super,2),'kd');text(center(:,1),center(:,2),num2str(mass_sp));
subplot(max_level,2,2*n_level+2)
scatter(X(:,1),X(:,2),5,map1D(super_potential,[lb,ub]));
        %%==========================================
        ind_track(n_level)=ind_super;
        if(~ismember(X(ind_mbp,:),Xsub,'rows'))
            disp('index outside of the targeting cluster')
            % break;
        end
        
    else
        disp('success')
        break;
    end
    n_level=n_level+1;
end
figure(findobj('Tag','ml'));subplot(nrow,1,1);hold on; h3=plot(X(ind_super,1),X(ind_super,2),'bo');
legend([h1,h2,h3],'global-MBP','first-super-MBP','last-super-MBP');




