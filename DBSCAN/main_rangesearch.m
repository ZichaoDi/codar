more off;
global epsm dimen N_cluster epsilon
dimen=3; %%==dimension of halo
files=dir('/homes/erangel/plank_halos_cf/*.bin');
for file=4;%11:20%1:length(files)
    %%===Load halo looping from Steve's database 
    [pathstr,name,ext]=fileparts(['/homes/erangel/plank_halos_cf/',files(file).name]);
    fid=fopen([pathstr,'/', name, ext]);
    data=fread(fid,'single');
    data=reshape(data,size(data,1)/3,3);
    N_particle=size(data,1); %%==total number of particles in the current halo
    fclose(fid);
    
    NN= N_particle;%[100 1000 2000 5000 10000 20000 30000 40000 50000];% 60000 70000 80000 N_particle] ;%%===downsample particles if NN < N_particle;
    count=[];
    err=zeros(length(NN),1);
    plotBPmap=0; %%== if 1 plot BP map on each level; 0 no plot.
    epsm=0.003; %%==distance threshold to determine if two particles should be combined
    t1=zeros(length(NN),1); %%==timing for brute force
    t3=zeros(length(NN),1); %%==timing for hierarchy
    rng('default')
    for i=1:length(NN)
        subN=NN(i)
        ind=unique(ceil(rand(subN,1)*N_particle));%%==uniformly downsample particles based on indices of each particle.
        NN(i)=length(ind);
        subN=NN(i);
        X=data(ind,1:dimen); %%==current full set of particles
        %%=========================brute force=================
        % tic;
        % [ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm); 
        % t1(i)=toc;
        % fprintf('time elapsed for global MBP: %d\n',t1(i));
        %%========================================================
        tic;
        tree = kd_buildtree(X,0); %%==build tree based on the full X;
        ind_super=0;
        epsilon=0.3;%%==initial distance threshold to rule out seeds too close to each other
        epsilon1=epsilon/2;%%==initial distance threshold for neigbhoring searching radius
%%%==================================== Global BP v.s. Local BP
mass=ones(subN,1);
potential=zeros(subN,1);
local_p=zeros(subN,1);
IDX=rangesearch(tree,X,epsilon1);%%==indices of particles as neighbors of the seeds
for i1=1:subN
    sub_ind=[1:i1-1,i1+1:subN];
    potential(i1)=-sum(mass(sub_ind)./max(sqrt(sum(bsxfun(@minus,X(i1,:),X(sub_ind,:)).^2,2)),epsm));
    local_p(i1)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(i1,:),X(IDX{i1},:)).^2,2)),epsm));%%==BP of each seed
end
return;
%%%================================================================================
        subind=[];
        count_t=1;
        for level=1:10 %%==recursive levels
            if(level==1)
                N_cluster=10; %%==initial number of seeds
                ind_k=unique(ceil(rand(N_cluster,1)*subN));
                N_cluster=length(ind_k);
                seeds_old=[];
                seeds=X(ind_k,:);%%==initial seeds
                super_potential=zeros(N_cluster,1);
                if N_cluster>1
                    IDX1=rangesearch(tree,seeds,epsilon1);%%==indices of particles as neighbors of the seeds
                    for j=1:N_cluster
                        super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,seeds(j,:),X(IDX1{j},:)).^2,2)),epsm));%%==BP of each seed
                        count{count_t}=[repmat(ind_k(j),length(IDX1{j}),1),IDX1{j}'];
            count_t=count_t+1;
                    end
                end
                sub2n=N_cluster;
                ind_super_diff=1;
                N_seeds_diff=N_cluster;
                clear ind_k IDX_temp;
            else
                N_cluster_old=N_cluster;
                N_cluster=min([N_cluster_old,length(smin)]);
                seeds_old=seeds;
                seeds=Xsub(smin(1:N_cluster),:); %%==assign current seeds as the BP peaks from previous level
                IDX1=rangesearch(tree,seeds,epsilon1);%%==indices of seeds' neighbors
                epsilon1=epsilon1/2; %%==recursively descrease the search radius
                totalind=[];
                super_potential_total=[];
                for j=1:N_cluster
                    subind=IDX1{j};%%==neighbors of current seed
                    sub2n=length(subind);
                    super_potential=zeros(sub2n,1);
                    count{count_t}=[];
                    for j1=1:sub2n
                        ind_other=setdiff(subind,subind(j1));
                        super_potential(j1)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(subind(j1),:),X(ind_other,:)).^2,2)),epsm));
                        count{count_t}=[count{count_t};repmat(subind(j1),length(ind_other),1), ind_other'];
                    end %%==BP for each particle in the neighborhood of current seed
                        count_t=count_t+1;
                    totalind=[totalind subind]; %%==combine considered particles
                    super_potential_total=[super_potential_total;super_potential];%%==combine considered paricles' BP
                end;
                subind=totalind;
                super_potential=super_potential_total;
                sub2n=length(subind);
                clear IDX_temp totalind super_potential_total;
                [~,sortind1]=sort(super_potential);
                ind_super_old=ind_super;
                ind_super=subind(sortind1(1)); %%==MBP provided by the current BP map
                ind_super_diff=norm(ind_super_old-ind_super);
                N_seeds_diff=N_cluster-N_cluster_old;
            end
            fprintf('Level %d, %d sampling points with %d particles\n',level,N_cluster,sub2n);
            %%===========plot BP map=========================================
            doplot(plotBPmap, X, real_potential,seeds,level,super_potential,subind);
            %%==============================================================
            [iconv,err1]=cnvtst(X,ind_super,ind_super_diff, seeds, seeds_old,ind_mbp,sub2n,N_seeds_diff);%%==test convergence; 
            if(iconv)
                err(i)=err1;
                break;
            end
            %%==Find peaks based on current BP map
            if(level==1)
                if(size(seeds,1)>1)
                    smin=FindPeak([seeds,super_potential],epsilon);
                else
                    smin=1;
                end
                Xsub=seeds;
            else
                smin=FindPeak([X(subind,:),super_potential],epsilon);
                Xsub=X(subind,:);
            end
            epsilon=epsilon/2;
        end
        t3(i)=toc;
        fprintf('time elapsed for finding super-MBP: %d\n',t3(i));
        doplot(plotBPmap,X,ind_mbp,ind_super,subN);
    end
end
