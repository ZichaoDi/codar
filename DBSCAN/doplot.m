function doplot(plotBPmap, X, varargin)
%%=========performance plot;
if(plotBPmap)
if(nargin==5)
    ind_mbp=varargin{1};
    ind_super=varargin{2};
    subN=varargin{3};
    figure,
    plot3(X(:,1),X(:,2),X(:,3),'r.','MarkerSize',0.1);hold on; 
    plot3(X(ind_mbp,1),X(ind_mbp,2),X(ind_mbp,3),'k+',X(ind_super,1),X(ind_super,2),X(ind_super,3),'bo','MarkerSize',15);
    legend('particle','global-mbp','super-mbp');title(num2str(subN));
    drawnow;
else
    real_potential=varargin{1};
    lb=min(varargin{1});
    ub=max(varargin{1});
    point_temp=varargin{2};
    level=varargin{3};
    super_potential=varargin{4};
    subind=varargin{5};
    figure,
    subplot(2,1,1)
    scatter3(X(:,1),X(:,2),X(:,3),1,map1D(real_potential,[lb,ub])); colorbar; hold on; 
    plot3(point_temp(:,1),point_temp(:,2),point_temp(:,3),'k+','MarkerSize',15)
    l=axis;
    title(['level ',num2str(level),', global BP'])
    subplot(2,1,2)
    if(level==1)
        scatter3(point_temp(:,1),point_temp(:,2),point_temp(:,3),1,map1D(super_potential,[lb,ub])); colorbar;
    else
        scatter3(X(subind,1),X(subind,2),X(subind,3),1,map1D(super_potential,[lb,ub])); colorbar;
    end
    axis(l);
    title('super BP')
    drawnow;
end
end
