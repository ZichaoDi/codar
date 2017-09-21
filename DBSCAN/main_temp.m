
            N_cluster_old=N_cluster;
            coreNo=5;
            N_cluster=min([coreNo,length(smin)]);
            point_temp_old=point_temp;
            point_temp=Xsub(smin(1:N_cluster),:);
            IDX1=rangesearch(tree,point_temp,epsilon);
            totalind=[];
            for j=1:N_cluster
                subind=IDX1{j};
                sub2n=length(subind);
                for j_sub=1:sub2n
                    ind_other{subind(j_sub)}=[ind_other{subind(j_sub)} setdiff(subind,subind(j_sub))];
                end
                totalind=[totalind subind];
            end;
            subind=unique(totalind);
            sub2n=length(subind);
            super_potential=zeros(sub2n,1);
            for j=1:sub2n
                super_potential(j)=-sum(1./max(sqrt(sum(bsxfun(@minus,X(totalind(j),:),X(unique(ind_other{j}),:)).^2,2)),epsm));
            end
