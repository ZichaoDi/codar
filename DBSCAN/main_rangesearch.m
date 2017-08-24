%% Load Data
more off;
global nrow epsm dimen
dimen=3;
rng('default')
files=dir('/homes/erangel/plank_halos_cf/*.bin');
figure,
for file=52%:length(files)
    [pathstr,name,ext]=fileparts(['/homes/erangel/plank_halos_cf/',files(file).name]);
    % N_particle=dir([pathstr,'/', name, ext]).bytes/3/4; % 3: x,y,z-axis; 4:single precision
    fid=fopen([pathstr,'/', name, ext]);
    data=fread(fid,'single');
    data=reshape(data,size(data,1)/3,3);
    N_particle=size(data,1);
    % plot3(data(:,1),data(:,2),data(:,3),'r.')
end
fclose(fid);
%%========================================================
NN=10000;%[100 1000 5000 10000 15000];% 20000 30000 40000 50000 60000 70000 80000 N_particle] ;% 
plotBPmap=1;
epsm=0.003;
t1=zeros(length(NN),1);
t3=zeros(length(NN),1);
for i=1:length(NN)
    i
    subN=NN(i);
    ind=unique(ceil(rand(subN,1)*N_particle));
    NN(i)=length(ind);
    subN=NN(i);
    X=data(ind,1:dimen);
    omega_global=[min(X(:,1:dimen));max(X(:,1:dimen))];
    tree = kd_buildtree(X,0);
    tic;
    [ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm);
    t1(i)=toc;
    fprintf('time elapsed for global MBP: %d\n',t1(i));
    lb=min(real_potential);
    ub=max(real_potential);
    tic;
    ind_super=0;
    epsb=0.2;
    epsilon=0.2;%sqrt(2)/2*epsb;
    subind=[];
    global_center=centroid(X);
    for level=1:10
        MinPts=0;
        if(level==1)
            N_cluster=150;
            ind_k=ceil(rand(N_cluster,1)*subN);
            point_temp=X(ind_k,:);
            % direction=[1 -1;1 1; -1 1; -1 -1; 1 0; 0 1; -1 0; 0 -1; 0 0];
            % % direction=direction(end,:);
            % point_temp=repmat(global_center,N_cluster,1)+epsb.*direction;
            super_potential=zeros(N_cluster,1);
            if N_cluster>1
                IDX1=rangesearch(tree,point_temp,epsilon);
                for j=1:N_cluster
                     ind_other=IDX1{j};
                     super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,point_temp(j,:),X(ind_other,:)).^2,2)),epsm));
                end
            end
            %%%==============================================
            % for j=1:N_cluster
            %     IDX_temp(IDX1{j})=j;
            % end
            % subind=find(IDX_temp~=0);
            % sub2n=length(subind);
            % super_potential=zeros(sub2n,1);
            % for j=1:sub2n
            %      ind_other=find(IDX_temp==IDX_temp(subind(j)));
            %      super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j),:),X(ind_other,:)).^2,2)),epsm));
            % end
            %%%==============================================
            clear ind_k IDX_temp;
            epsilon1=epsilon;
            fprintf('Level %d, %d sampling points\n',level,N_cluster);
        else
            N_cluster_old=N_cluster;
            N_cluster=min([4,length(smin)]);
            point_temp_old=point_temp;
            point_temp=[xq(smin(1:N_cluster)),yq(smin(1:N_cluster)),zq(smin(1:N_cluster))];
            IDX1=rangesearch(tree,point_temp,epsilon1);
            % super_potential_inter=zeros(N_cluster,1);
            % epsilon1=epsilon1/2;
            totalind=[];
            super_potential_total=[];
            for j=1:N_cluster
                % IDX1=rangesearch(tree,point_temp(j,:),epsilon1);
                % IDX_temp{j}=IDX1{1};
                % ind_other=IDX1{1};
                % super_potential_inter(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,point_temp(j,:),X(ind_other,:)).^2,2)),epsm));
                subind=IDX1{j};
                sub2n=length(subind);
                super_potential=zeros(sub2n,1);
                for j=1:sub2n
                    ind_other=setdiff(subind,j);
                    super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j),:),X(ind_other,:)).^2,2)),epsm));
                end
                totalind=[totalind subind];
                super_potential_total=[super_potential_total;super_potential];
            end;
            super_potential=super_potential_total;
            subind=totalind;
            % [~,ind_super_inter]=min(super_potential_inter);
            % subind=IDX_temp{ind_super_inter};%find(IDX(ind_super_inter,:)==1);
            clear IDX_temp;
            % sub2n=length(subind);
            % super_potential=zeros(sub2n,1);
            % for j=1:sub2n
            %     ind_other=setdiff(subind,j);
            %     super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j),:),X(ind_other,:)).^2,2)),epsm));
            % end
            fprintf('Level %d, %d sampling points with %d particles\n',level,N_cluster,sub2n);
            [~,sortind1]=sort(super_potential);
            ind_super_old=ind_super;
            if(length(subind)>0)
            ind_super=subind(sortind1(1));
            if(norm(ind_super_old-ind_super)==0 | norm(point_temp_old(1:N_cluster,:)-point_temp)<1e-3)
                if(N_cluster==N_cluster_old)
                    if(norm(point_temp_old(1:N_cluster,:)-point_temp)<1e-3)
                        disp('no difference for seeds')
                    end
                end
                if(ind_super==ind_mbp)
                    disp('converge to global MBP')
                else
                    err=norm(X(ind_mbp,:)-X(ind_super,:));
                    fprintf('converge with error %e\n',err);
                end
                break;
            end
        else
            fprintf('no points for the sub region\n');
            break;
        end

        end
        %%===========================================================================
        if(plotBPmap)
            figure,
            subplot(2,1,1)
            scatter3(X(:,1),X(:,2),X(:,3),5,map1D(real_potential,[lb,ub])); colorbar;
            hold on; plot3(point_temp(:,1),point_temp(:,2),point_temp(:,3),'k+')
            l=axis;
            title(['level ',num2str(level),', global BP'])
            subplot(2,1,2)
            if(level==1)
                scatter3(point_temp(:,1),point_temp(:,2),point_temp(:,3),5,map1D(super_potential,[lb,ub])); colorbar;
            else
                scatter3(X(subind,1),X(subind,2),X(subind,3),5,map1D(super_potential,[lb,ub])); colorbar;
            end
            axis(l);
            title('super BP')
        end
            %%===========================================================================
        if(level==1)
            if(size(point_temp,1)>1)
                omega=[min(point_temp(:,1:dimen)); max(point_temp(:,1:dimen))];
                gridsize=10;
                F=scatteredInterpolant(point_temp(:,1),point_temp(:,2),point_temp(:,3),super_potential);
                [xq,yq,zq]=meshgrid(linspace(omega(1,1),omega(2,1),gridsize),linspace(omega(1,2),omega(2,2),gridsize),linspace(omega(1,3),omega(2,3),gridsize));
                % vq=griddata(point_temp(:,1),point_temp(:,2),point_temp(:,3),super_potential,xq,yq,zq);
                vq=F(xq,yq,zq);%interp3(point_temp(:,1),point_temp(:,2),point_temp(:,3),super_potential,xq,yq,zq);
                vq(isnan(vq))=0;
                BW=imregionalmax(-vq);
                smin=find(BW(:)==1);
            else
                smin=1;
                xq=point_temp(1);yq=point_temp(2);
            end
        else
            omega=[min(X(subind,1:dimen)); max(X(subind,1:dimen))];
            gridsize=10;
            [xq,yq,zq]=meshgrid(linspace(omega(1,1),omega(2,1),gridsize),linspace(omega(1,2),omega(2,2),gridsize),linspace(omega(1,3),omega(2,3),gridsize));
            F=scatteredInterpolant(X(subind,1),X(subind,2),X(subind,3),super_potential);
            vq=F(xq,yq,zq);
            vq(isnan(vq))=0;
            BW=imregionalmax(-vq);
            smin=find(BW(:)==1);
        end
    end
    ind_super=subind(sortind1(1:min(1,length(sortind1))));
    t3(i)=toc;
    % [~,sortind]=sort(real_potential);
    % in1=[]; for j=1:sub2n,in1(j)=find(sortind==subind(j));end
    fprintf('time elapsed for finding super-MBP: %d\n',t3(i));
    figure,plot3(X(:,1),X(:,2),X(:,3),'r.');hold on; plot3(X(ind_mbp,1),X(ind_mbp,2),X(ind_mbp,3),'k+',X(ind_super,1),X(ind_super,2),X(ind_super,3),'bo',global_center(1),global_center(2),global_center(3),'gd');legend('particle','global-mbp','super-mbp','centroid');title(num2str(subN));
    
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
