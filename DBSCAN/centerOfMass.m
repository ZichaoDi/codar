function centroid = centerOfMass(A,varargin)
centroid=sum(A,1)/size(A,1);

% %% INPUT CHECK
% narginchk(0,1);
% nargoutchk(0,1);
% fname = 'centerOfMass';
% 
% % Checked required inputs
% validateattributes(A,{'numeric'},{'real','finite'},fname,'A',1);
% 
% %% INITIALIZE VARIABLES
% A(isnan(A)) = 0;
% if ~(strcmpi(class(A),'double') || strcmpi(class(A),'single'))
%     A = single(A);
% end
% if any(A(:)<0)
%     warning('MATLAB:centerOfMass:neg','Array A contains negative values.');
% end
% 
% %% PROCESS
% sz = size(A);
% nd = ndims(A);
% M = sum(A(:));
% C = zeros(1,nd);
% if M==0
%     C = [];
% else
%     for ii = 1:nd
%         shp = ones(1,nd);
%         shp(ii) = sz(ii);
%         rep = sz;
%         rep(ii) = 1;
%         ind = repmat(reshape(1:sz(ii),shp),rep);
%         C(ii) = sum(ind(:).*A(:))./M;
%     end
% end
% 
% % Assemble the VARARGOUT cell array
% varargout = {C};
% 
% end % MAIN
