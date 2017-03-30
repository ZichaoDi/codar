function [ind,potential]=mbp(X,mass);
N=size(X,1);
potential=zeros(N,1);
dist=ones(N,N);
for i=1:N
    for j=1:N
        if(i~=j)
            dist=norm(X(i,:)-X(j,:));
            potential(i)=potential(i)-mass(j)/(dist+0.3);
            % if(dist>=0.05)
            %     potential(i)=potential(i)-mass(j)/dist;
            % end
        end
    end
end
[~,ind]=min(potential);

