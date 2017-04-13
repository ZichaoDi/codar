%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function h1=PlotClusterinResult(X, IDX, ind_mbp, epsilon, MinPts,clustering)
global epsm
    figure, 
    subplot(2,1,1);plot(X(:,1),X(:,2),'r.');
    hold on; h1=plot(X(ind_mbp,1),X(ind_mbp,2),'k*');
    title(num2str(epsm))
    legend(h1,'global-MBP');

    axis equal;
    subplot(2,1,2);

    k=max(IDX);

    Colors=hsv(k);

    Legends = {};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = '.';
            MarkerSize = 8;
            Color = Colors(i,:);
            Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 6;
            Color = [0 0 0];
            if ~isempty(Xi)
                Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);
        end
        hold on;
    end
    hold off;
    axis equal;
    grid on;
    legend(Legends);
    legend('Location', 'NorthEastOutside');
    if(strcmp(clustering,'dbscan'))
        title(['DBSCAN (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
    elseif(strcmp(clustering,'kmeans'))
        title('Kmeans');
    end

end
