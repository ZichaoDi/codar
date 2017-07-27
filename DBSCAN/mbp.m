function [ind,potential]=mbp(X,mass,epsm);
N=size(X,1);
potential=zeros(N,1);
for i=1:N
    sub_ind=[1:i-1,i+1:N];
    potential(i)=-sum(mass(sub_ind)./max(sqrt(sum(bsxfun(@minus,X(i,:),X(sub_ind,:)).^2,2)),epsm));
end
% D = bsxfun(@plus,dot(X,X,1),dot(X,X,1)')-2*(X'*X);
[~,ind]=min(potential);

