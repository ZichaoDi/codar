function IDX = rangesearch(tree,point,epsilon);
global dimen
range=[-epsilon*ones(1,dimen); epsilon*ones(1,dimen)];
IDX=cell(size(point,1),1);
for k=1:size(point,1)
IDX_temp  = kd_rangequery(tree,point(k,:)',range);
IDX{k}=IDX_temp';
end

