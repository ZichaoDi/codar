%% Load Data
more off;
global nrow epsm dimen coreNo epsilon
dimen=3;
files=dir('/homes/erangel/plank_halos_cf/*.bin');
% figure,
starting=101;
for file=starting:length(files)
    [pathstr,name,ext]=fileparts(['/homes/erangel/plank_halos_cf/',files(file).name]);
    fid=fopen([pathstr,'/', name, ext]);
    data=fread(fid,'single');
    data=reshape(data,size(data,1)/3,3);
    N_particle=size(data,1);
    % plot3(data(:,1),data(:,2),data(:,3),'r.')
fclose(fid);
%%========================================================
NN= N_particle;% 
err=zeros(length(NN),1);
plotBPmap=0;
epsm=0.003;
t1=zeros(length(NN),1);
t3=zeros(length(NN),1);
rng('default')
for i=1:length(NN)
    i
    subN=NN(i);
    ind=1:N_particle;%randperm(subN);%unique(ceil(rand(subN,1)*N_particle));
    % ind=randperm(subN);%unique(ceil(rand(subN,1)*N_particle));
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
    epsilon=0.2;%sqrt(2)/2*epsb;
    subind=[];
    global_center=centroid(X);
    for level=1:10
        MinPts=0;
        if(level==1)
            N_cluster=50;
            ind_k=unique(ceil(rand(N_cluster,1)*subN));
            N_cluster=length(ind_k);
            point_temp=X(ind_k,:);
            super_potential=zeros(N_cluster,1);
            if N_cluster>1
                IDX1=rangesearch(tree,point_temp,epsilon);
                for j=1:N_cluster
                     ind_other=IDX1{j};
                     super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,point_temp(j,:),X(ind_other,:)).^2,2)),epsm));
                end
            end
            clear ind_k IDX_temp;
            fprintf('Level %d, %d sampling points\n',level,N_cluster);
        else
            N_cluster_old=N_cluster;
            coreNo=5;
            N_cluster=min([coreNo,length(smin)]);
            point_temp_old=point_temp;
            % point_temp=[xq(smin(1:N_cluster)),yq(smin(1:N_cluster)),zq(smin(1:N_cluster))];
            point_temp=Xsub(smin(1:N_cluster),:);
            IDX1=rangesearch(tree,point_temp,epsilon);
            totalind=[];
            super_potential_total=[];
            for j=1:N_cluster
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
            subind=totalind;
            super_potential=super_potential_total;
            sub2n=length(subind);
            clear IDX_temp;
            fprintf('Level %d, %d sampling points with %d particles\n',level,N_cluster,sub2n);
            [~,sortind1]=sort(super_potential);
            ind_super_old=ind_super;
        if(sub2n>1)
            ind_super=subind(sortind1(1));
            if(norm(ind_super_old-ind_super)==0)
                if(N_cluster==N_cluster_old)
                    if(norm(point_temp_old(1:N_cluster,:)-point_temp)<1e-3)
                        disp('no difference for seeds')
                    end
                end
                if(ind_super==ind_mbp)
                    err(i)=0;
                    disp('converge to global MBP')
                else
                    err(i)=norm(X(ind_mbp,:)-X(ind_super,:));
                    fprintf('converge with error %e\n',err(i));
                end
                break;
            end
        elseif(sub2n==1)
            ind_super=subind;
            err(i)=norm(X(ind_mbp,:)-X(ind_super,:));
            fprintf('converge with error %e\n',err(i));
            break;
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
        gridsize=20;
        if(level==1)
            if(size(point_temp,1)>1)
                % omega=[min(point_temp(:,1:dimen)); max(point_temp(:,1:dimen))];
                % F=scatteredInterpolant(point_temp(:,1),point_temp(:,2),point_temp(:,3),super_potential);
                % [xq,yq,zq]=meshgrid(linspace(omega(1,1),omega(2,1),gridsize),linspace(omega(1,2),omega(2,2),gridsize),linspace(omega(1,3),omega(2,3),gridsize));
                % % vq=F(xq,yq,zq);
                % vq=griddata(point_temp(:,1),point_temp(:,2),point_temp(:,3),super_potential,xq,yq,zq);
                % vq(isnan(vq))=0;
                % BW=imregionalmax(-vq);
                % smin=find(BW(:)==1);
                smin=FindPeak([point_temp,super_potential]);
            else
                smin=1;
                xq=point_temp(1);yq=point_temp(2);
            end
            Xsub=point_temp;

        else
            % omega=[min(X(subind,1:dimen)); max(X(subind,1:dimen))];
            % [xq,yq,zq]=meshgrid(linspace(omega(1,1),omega(2,1),gridsize),linspace(omega(1,2),omega(2,2),gridsize),linspace(omega(1,3),omega(2,3),gridsize));
            % F=scatteredInterpolant(X(subind,1),X(subind,2),X(subind,3),super_potential);
            % vq=F(xq,yq,zq);
            % vq(isnan(vq))=0;
            % BW=imregionalmax(-vq);
            % smin=find(BW(:)==1);
            smin=FindPeak([X(subind,:),super_potential]);
            Xsub=X(subind,:);
        end
    end
    ind_super=subind(sortind1(1:min(1,length(sortind1))));
    t3(i)=toc;
    fprintf('time elapsed for finding super-MBP: %d\n',t3(i));
    if(plotBPmap)
        figure,plot3(X(:,1),X(:,2),X(:,3),'r.');hold on; plot3(X(ind_mbp,1),X(ind_mbp,2),X(ind_mbp,3),'k+',X(ind_super,1),X(ind_super,2),X(ind_super,3),'bo',global_center(1),global_center(2),global_center(3),'gd');legend('particle','global-mbp','super-mbp','centroid');title(num2str(subN));
    end
end
errHalo(file-starting+1)=err(end);
tg(file-starting+1)=t1;
ts(file-starting+1)=t3;
N_tot(file-starting+1)=N_particle;
save(['realHalo',num2str(starting),'.mat'],'N_tot','tg','ts','errHalo');

end
