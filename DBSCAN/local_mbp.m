global epsm dimen N_cluster epsilon
dimen=3; %%==dimension of halo
files=dir('/homes/erangel/plank_halos_cf/*.bin');
for file=1:length(files)
    %%===Load halo looping from Steve's database 
    [pathstr,name,ext]=fileparts(['/homes/erangel/plank_halos_cf/',files(file).name]);
    fid=fopen([pathstr,'/', name, ext]);
    data=fread(fid,'single');
    data=reshape(data,size(data,1)/3,3);
    N_particle=size(data,1); %%==total number of particles in the current halo
    fclose(fid);
    
    NN= [100 1000 2000 5000 10000 20000 30000 40000 50000 60000 70000 80000 N_particle] ;%%===downsample particles if NN < N_particle;
    epsm=0.003; %%==distance threshold to determine if two particles should be combined
    rng('default')
    for i=1:length(NN)
        subN=NN(i)
        ind=unique(ceil(rand(subN,1)*N_particle));%%==uniformly downsample particles based on indices of each particle.
        NN(i)=length(ind);
        subN=NN(i);
        X=data(ind,1:dimen); %%==current full set of particles
        %%=========================brute force=================
        % tic;
        [ind_mbp,real_potential]=mbp(X,ones(subN,1),epsm); 
        mbp_file{(file-1)*length(NN)+i}=[ind;ind_mbp];
        save mbp_file mbp_file
        % t1(i)=toc;
        % fprintf('time elapsed for global MBP: %d\n',t1(i));
        % %%========================================================
        % tree = kd_buildtree(X,0); %%==build tree based on the full X;
        % eps1=linspace(0,1,100);
        % err=zeros(length(eps1),1);
        % for eps_i=1:length(eps1)
        %     tic;
        %     IDX1=rangesearch(tree,X,eps1(eps_i));%%==indices of particles as neighbors of the seeds
        %     for j=1:subN
        %         super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(j,:),X(IDX1{j},:)).^2,2)),epsm));%%==BP of each seed
        %     end
        %     [~,ind_super]=min(super_potential);
        %     err(eps_i)=norm(X(ind_mbp,:)-X(ind_super,:));
        %     t2(eps_i)=toc;
        % end
    end
end
% save(['local_mbp',num2str(NN),'.mat'],'t2','t1','err');
