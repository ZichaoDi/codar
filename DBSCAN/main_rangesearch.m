%% Load Data
more off;
global nrow epsm
rng('default')
[pathstr,name,ext]=fileparts('/homes/erangel/plank_halos_cf/499455452.bin');
% N_particle=dir([pathstr,'/', name, ext]).bytes/3/4; % 3: x,y,z-axis; 4:single precision
fid=fopen([pathstr,'/', name, ext]);
data=fread(fid,'single');
data=reshape(data,size(data,1)/3,3);
N_particle=size(data,1);
% figure,plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
% %%========================================================
% subN=10000;%size(data,1);
% ind=ceil(rand(subN,1)*N_particle);
% X=data(ind,1:2);
% %%========================================================
NN=[10 100 1000 5000 10000 15000 ];%20000 30000 40000 50000 60000 70000 80000 N_particle] ;
epsm=0.003;
t1=zeros(length(NN),1);
t3=zeros(length(NN),1);
for i=1:length(NN)
    i
    subN=NN(i);
    ind=ceil(rand(subN,1)*N_particle);
    X=data(ind,1:2);
    global_center=centroid(X);
    tic;
    [ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm);
    t1(i)=toc;
    fprintf('time elapsed for global MBP: %d\n',t1(i));
    lb=min(real_potential);
    ub=max(real_potential);
    tic;
    ind_super=0;
    epsilon=0.3;
    IDX=[];
    subind=[];
    ind_k=[];
    for level=1:10
        MinPts=0;
        N_cluster=5;
        if(level==1)
            ind_k=ceil(rand(N_cluster,1)*subN);
            point_temp=X(ind_k,:);
            IDX1=rangesearch(X,X(ind_k,:),epsilon);
            for j=1:N_cluster
                IDX(IDX1{j})=j;
            end
            subind=find(IDX~=0);
            sub2n=length(subind);
            super_potential=zeros(sub2n,1);
            for j=1:sub2n
                 ind_other=find(IDX==IDX(subind(j)));
                 super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j),:),X(ind_other,:)).^2,2)),epsm));
            end
        else
            N_cluster=min([sub2n,5,length(smin)]);
            ind_k=sortind1(1:N_cluster);
            N_cluster=length(ind_k);
            point_temp=[xq(smin(1:N_cluster)),yq(smin(1:N_cluster))];
            IDX=sparse(N_cluster,subN);
            super_potential_inter=zeros(N_cluster,1);
            for j=1:N_cluster
                IDX1=rangesearch(X,point_temp(j,:),epsilon);
                IDX(j,IDX1{1})=1;
                ind_other=IDX1{1};
                super_potential_inter(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,point_temp(j,:),X(ind_other,:)).^2,2)),epsm));
            end;
            [~,ind_super_inter]=min(super_potential_inter);
            subind=find(IDX(ind_super_inter,:)==1);
            sub2n=length(subind);
            super_potential=zeros(sub2n,1);
            for j=1:sub2n
                ind_other=setdiff(subind,j);
                super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j),:),X(ind_other,:)).^2,2)),epsm));
            end
        end
        fprintf('Level %d, # of sub-clusters: %d \n',level,N_cluster);
        [~,sortind1]=sort(super_potential);
        ind_super_old=ind_super;
        ind_super=subind(sortind1(1));
        %%===========================================================================
        % figure,
        % subplot(2,1,1)
        % scatter(X(:,1),X(:,2),5,map1D(real_potential,[lb,ub])); colorbar;
        % hold on; plot(point_temp(:,1),point_temp(:,2),'k+')
        % l=axis;
        % title('global BP')
        % subplot(2,1,2)
        % scatter(X(subind,1),X(subind,2),5,map1D(super_potential,[lb,ub])); colorbar;
        % axis(l);
        % title('super BP')
        %%===========================================================================
        if(norm(ind_super_old-ind_super)==0)
            disp('converge')
            break;
        end
        omega=[min(X(subind,1)) max(X(subind,1)) min(X(subind,2)) max(X(subind,2))];
        [xq,yq]=meshgrid(linspace(omega(1),omega(2),100),linspace(omega(3),omega(4),100));
        vq=griddata(X(subind,1),X(subind,2),super_potential,xq,yq);
        [xymax,smax,xymin,smin]=extrema2(vq);
    end
    ind_super=subind(sortind1(1:min(1,length(sortind1))));
    t3(i)=toc;
    [~,sortind]=sort(real_potential);
    % in1=[]; for j=1:sub2n,in1(j)=find(sortind==subind(j));end
    fprintf('time elapsed for finding super-MBP: %d\n',t3(i));
    figure,plot(X(:,1),X(:,2),'r.');hold on; plot(X(ind_mbp,1),X(ind_mbp,2),'k+',X(ind_super,1),X(ind_super,2),'bo',global_center(1),global_center(2),'gd');legend('particle','global-mbp','super-mbp','centroid');title(num2str(subN));
    
    % nrow=2;
    % h1=PlotClusterinResult(X,IDX,ind_mbp,epsilon,MinPts,clustering);
    % figure(findobj('Tag','cluster'));
    % hold on; plot(X(ind_mbp,1),X(ind_mbp,2),'k+',X(ind_super,1),X(ind_super,2),'kd');
    % figure(findobj('Tag','ml'));
    % subplot(nrow,1,1);hold on; h2=plot(X(ind_super,1),X(ind_super,2),'go');
    %
    % legend([h1,h2],'global-MBP','Super-MBP');
    % subplot(nrow,1,nrow), plot(real_potential(sortind),'r.-');hold on; plot(in1,map1D(super_potential,[lb,ub]),'b.');hold off;
    % legend('real-potential','super-potential')
    % % subplot(nrow,1,4), plot(real_potential(sortind1),'r.-');hold on; plot(map1D(super_potential(sortind1),[lb,ub]),'b.-');hold off;
    % legend('real-potential','super-potential')
end
