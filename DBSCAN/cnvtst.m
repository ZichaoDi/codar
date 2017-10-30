function [iconv,err1]=cnvtst(X,ind_super,ind_super_diff, point_temp, point_temp_old,ind_mbp,sub2n,N_cluster_diff)
%%==convergence test
iconv=0;
err1=[];
if(sub2n>1)
    if(ind_super_diff==0)
        if(N_cluster_diff==0)
            if(norm(point_temp_old-point_temp)<1e-3)
                disp('no difference for seeds')
            end
        end
        if(ind_super==ind_mbp)
            err1=0;
            disp('converge to global MBP')
        else
            err1=norm(X(ind_mbp,:)-X(ind_super,:));
            fprintf('converge with error %e\n',err1);
        end
        iconv=1;
    end
elseif(sub2n==1)
    err1=norm(X(ind_mbp,:)-X(ind_super,:));
    fprintf('converge with error %e\n',err1);
    iconv=1;
else
    fprintf('no points for the sub region\n');
    iconv=1;
end
