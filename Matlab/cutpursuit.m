function in_component=cutpursuit(pts,K,regStrength,edgeStength,mode,speed,verbose,edgeWeight,Euv,nodeWeight)
%CUTPURSUIT Segment 2D/3D/ND array based on L0 cut-pursuit algorithm
%   Copyright (c) 2018 Loic Landrieu
%   Reference: Cut Pursuit: fast algorithms to learn piecewise constant functions on general weighted graphs, L. Landrieu and G. Obozinski, SIAM Journal on Imaging Science 2017, Vol. 10, No. 4 : pp. 1724-1766 [hal link]
%   https://github.com/loicland/cut-pursuit

nNodes=size(pts,1);
if ~exist('Euv','var')
    Mdl = KDTreeSearcher(pts);
    [minIdxs,~] = knnsearch(Mdl,pts,'K',K);

    %graph
    Eu = (1:nNodes)';
    Euv=[];
    for j=1:K-1
        cEuv=[Eu,minIdxs(:,j+1)];
        Euv=[Euv;cEuv];
    end
    Euv=Euv-1;
end

if ~exist('edge_weight','var')
    edgeWeight = edgeStength*ones(size(Euv,1),1);
end

if ~exist('node_weight','var')
    nodeWeight = ones(nNodes,1);
end

[~, in_component, ~] = L0_cut_pursuit_segmentation(single(pts'), uint32(Euv(:,1)'), uint32(Euv(:,2)'), single(regStrength), single(edgeWeight), single(nodeWeight), mode, speed, verbose);
end
