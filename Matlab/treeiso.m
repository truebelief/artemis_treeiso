warning off
fclose all;
close all;
fprintf('Individual-tree isolation (treeiso) from terrestrial laser scanning point clouds: Batch mode\n');
fprintf('Ground points should be already removed\n');
fprintf('The University of Lethbridge - Artemis Lab\n');
fprintf('Copyright - Zhouxin Xi (zhouxin.xi@uleth.ca)\n');

%If you use your own TLS data, before running this code, it is recommended to clean the point clouds and remove the noise, e.g. using CloudCompare(Sober filter,NN=10,std=1)
%It is better to decimate the point clouds to ~2cm resolution before the code

fpath=input('Please enter path including all scan files (file name format: *.laz): ','s'); %(file name must start with a unique ID number)
foutpath='output';

if exist(foutpath,'dir')<1
    mkdir(foutpath);
end

PR_REG_STRENGTH1=1.0;%lambda1
PR_MIN_NN1=5;%K1:key parameter
PR_EDGE_STRENGTH1=1;

PR_REG_STRENGTH2=20;%lambda2:key parameter
PR_MIN_NN2=20;%K2:key parameter
PR_EDGE_STRENGTH2=1;

PR_REL_HEIGHT_LENGTH_RATIO=0.5;%rho
PR_VERTICAL_WEIGHT=0.5;%w:key parameter


    
PR_DECIMATE_RES1=0.05;
PR_DECIMATE_RES2=0.1;
PR_MAX_GAP=2.0;

option=1;
speed=0;
verbose=0;

fls=dir([fpath,'/*.laz']);
fnames={fls.name}';

for t=1:numel(fnames)
    [~,fname,~]=fileparts(fnames{t});
    fprintf('Processing (%d/%d): %s @%s',t,numel(fnames),fname,datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
    
    %% 1. initial 3D segmentation
    [hdr,s]  = las2mat(['-i ',strcat(fpath,'\',fname,'.laz')]);
  
    scenePts=[s.x,s.y,s.z,double(s.intensity)];
    segs=[scenePts(:,1:end),zeros(size(scenePts,1),1)];
        
    subPos=scenePts(:,1:3);
    if size(subPos,1)<PR_MIN_NN1*2
        continue;
    end
    subPos=subPos-mean(subPos,1);
    
    [~,ia,ic]=unique(floor((subPos-min(subPos,[],1))/PR_DECIMATE_RES1)+1,'rows');
    decPos=subPos(ia,:);

    nNodes=size(decPos,1);
    if nNodes>0
        inComponent=cutpursuit(decPos,PR_MIN_NN1,PR_REG_STRENGTH1,PR_EDGE_STRENGTH1,option,speed,verbose);
        segs(:,end)=double(inComponent(ic))+1;
    end

    %% 2. Bottom-up 2D segmentation
    segs=[segs,zeros(size(segs,1),1)];
    segsMin=min(segs(:,1:3),[],1);
    segs(:,1:3)=segs(:,1:3)-segsMin;
    
    clusterIdx=segs(:,end-1);
    [clusterU,clusterUIdx,clusterV]=unique(clusterIdx);
    [clusterVSorted,clusterVIdx]=sort(clusterV);
    clusterVGroup = mat2cell(clusterVIdx, diff([0;find([(diff(clusterVSorted) ~= 0);1])]),1);
 
    currentClusterDecs=zeros(size(segs,1),4)-1;
    currentClusterIs={};
    
    % extract each cluster centroid and decimate cluster into low-res points
    clusterCentroids=zeros(numel(clusterU),3);
    iter=0;
    for i=1:numel(clusterU)
        currentClusterPos=segs(clusterVGroup{i},1:3);
        [~,ia,ic]=unique(floor((currentClusterPos-min(currentClusterPos,[],1))/PR_DECIMATE_RES2)+1,'rows');

        startI=iter+1;
        iter=iter+numel(ia);
        end_i=iter;
        currentClusterDecs(startI:end_i,:)=[currentClusterPos(ia,:),repmat(i,numel(ia),1)];
        currentClusterIs{i}=startI:end_i;
        clusterCentroids(i,:)=mean(currentClusterPos,1);
    end
    currentClusterDecs(currentClusterDecs(:,1)==-1,:)=[];
    currentClusterDecsMean=mean(currentClusterDecs(:,1:3),1);
    currentClusterDecs(:,1:3)=currentClusterDecs(:,1:3)-currentClusterDecsMean;

    % for each centroid, extract its nn cluster point distance to be the edge weight
    Mdl = KDTreeSearcher(clusterCentroids);
    [minIdxsC,~] = knnsearch(Mdl,clusterCentroids,'K',PR_MIN_NN2);
    
    nnDists=zeros(size(minIdxsC))+10;
    for i=1:size(minIdxsC,1)
%         current_cluster=segs(clusterVGroup{minIdxsC(i,1)},1:3);
        currentClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,1)},1:3);
        
        for j=2:size(minIdxsC,2)
            nnClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,j)},1:3);
            [~,D]=jitknnsearch(nnClusterDec,currentClusterDec);
            nnDists(i,j)=min(D);
        end
    end
    nnDists(:,1)=0;
    

    % prepare pair-wise edge weight and segment (decimated points as graph nodes)
    nNodes=size(currentClusterDecs,1);%per point, not per cluster centroid
    Mdl = KDTreeSearcher(currentClusterDecs(:,1:3));
    [minIdxs,Ds] = knnsearch(Mdl,currentClusterDecs(:,1:3),'K',PR_MIN_NN2);
    
    edgeWeight=zeros(size(minIdxs,1)*size(minIdxs,2),1)-1;
    Euv=zeros(size(minIdxs,1)*size(minIdxs,2),2);
    for i=1:size(minIdxs,1)
        currentNode=currentClusterDecs(i,4);
        currentDists=nnDists(currentNode,:);
        for j=2:size(minIdxs,2)
            nnNode=currentClusterDecs(minIdxs(i,j),4);
            nnDist=currentDists(minIdxsC(currentNode,:)==nnNode);
            if numel(nnDist)>1
                nnDist=nnDist(1);
            end
            if nnDist<PR_MAX_GAP
                Euv((i-1)*size(minIdxs,2)+j,:)=[i,minIdxs(i,j)];
                edgeWeight((i-1)*size(minIdxs,2)+j)=10./(nnDist/0.01);
            end
        end
    end
    Euv(edgeWeight==-1,:)=[];
    edgeWeight(edgeWeight==-1)=[];
    
    Euv=Euv-1;
    inComponent=cutpursuit(currentClusterDecs(:,1:2),PR_MIN_NN2,PR_REG_STRENGTH2,PR_EDGE_STRENGTH2,option,speed,verbose,edgeWeight,Euv);
    for i=1:numel(clusterU)        
        currentComponent=mode(inComponent(currentClusterIs{i})+1);
        segs(clusterVGroup{i},end)=currentComponent;
    end


    %% 3. Global edge refinement
    % there will still be over-segmented crowns, not merged into one tree yet;
    segs=[segs,segs(:,end)];
    

    Euvw=[Euv+1,0.01./(edgeWeight./10.0)];
    clusterGroupMap=[clusterCentroids,segs(clusterUIdx,end-2:end-1)];
    
    [clusterMapU,~,clusterMapV] = unique(clusterGroupMap(:,end));
    [clusterMapVSorted,clusterMapVIdx]=sort(clusterMapV);
    clusterMapVGroup = mat2cell(clusterMapVIdx, diff([0;find([(diff(clusterMapVSorted) ~= 0);1])]),1);

    iter=1;
    
    %Stem identifier: ratio of elevation difference from neighbors to segment length < PR_REL_HEIGHT_LENGTH_RATIO
    toMergeIds=[0,0];
    prevToMergeIds=[0];
    mergedRemainIds=[];
    while (numel(toMergeIds)~=numel(prevToMergeIds) && numel(toMergeIds)>0)
        prevToMergeIds=toMergeIds;

        [groupU,~,groupV] = unique(segs(:,end));
        [groupVSorted,groupVIdx]=sort(groupV);
        groupVGroup = mat2cell(groupVIdx, diff([0;find([(diff(groupVSorted) ~= 0);1])]),1);
        nv = cellfun(@length, groupVGroup);
    
        nGroups=numel(groupVGroup);
        groupFeatures=zeros(nGroups,5);
        groupHulls=cell(nGroups,1);
        
        for i=1:numel(groupVGroup)
            groupPts=segs(groupVGroup{i},:);
            groupFeatures(i,1:3)=trimmean(groupPts(:,1:3),0.2);%centroid
            groupFeatures(i,4)=min(groupPts(:,3));%elevation
            groupFeatures(i,5)=max(groupPts(:,3))-min(groupPts(:,3));%length
            if size(groupPts,1)>3
                [groupHull,~] = convhull(groupPts(:,1),groupPts(:,2));
                groupHulls{i}=groupPts(groupHull,1:2);%convec hull to calculate overlap
            end
        end

        PR_MIN_NN3=20;
        nGroups=size(groupFeatures,1);%per point, not per cluster centroid
        Mdl = KDTreeSearcher(groupFeatures(:,1:2));
        [groupNNIdxC,groupNNCDs] = knnsearch(Mdl,groupFeatures(:,1:2),'K',PR_MIN_NN3);

        %mean centroid distance between segments
        sigmaD=mean(groupNNCDs(:,2));
                
        %if the segment elevation compared to neighboring segment elevations is high, and the segment is short, it is very likely to the be topper unmerged segments
        toMergeIds=zeros(nGroups,1);

        for i=1:nGroups        
                currentGroupFt=groupFeatures(i,4:5);
                nnGroupId=groupNNIdxC(i,1:end);
                nnGroupFt=groupFeatures(nnGroupId,4:5);
                currentGroupRelHt=(currentGroupFt(1)-min(nnGroupFt(:,1)))/currentGroupFt(2);
                if abs(currentGroupRelHt)>PR_REL_HEIGHT_LENGTH_RATIO
                    if (currentGroupFt(2)/median(nnGroupFt(:,2)))>1.5%set flags, for the first iteration, just merge longer segments first, for the remaining iterations, merge others
                        toMergeIds(i)=i;
                    else
                        toMergeIds(i)=-i;
                    end
                end
        end

        toMergeIds(toMergeIds==0)=[];
        remainIds=setdiff((1:nGroups)',abs(toMergeIds));
        if (iter==1) && sum(toMergeIds>=0)>0
            toMergeIds(toMergeIds<0)=[];
        else
            toMergeIds=abs(toMergeIds);
        end

        [groupNNIdx,~] = jitknnsearch(groupFeatures(toMergeIds,1:2),groupFeatures(remainIds,1:2),min(PR_MIN_NN3,numel(remainIds)));

        nToMergeIds=numel(toMergeIds);
        nRemainId=numel(remainIds);
        for i=1:nToMergeIds
            toMergeId=toMergeIds(i);
            currentClusterCentroids=clusterGroupMap(vertcat(clusterMapVGroup{toMergeId}),:);

            nNNs=size(groupNNIdx,2);
            filterMetrics=zeros(nNNs,5);%horizonal 2D overlapping ratio, vertical 1D overlapping ratio, 3D NN cluster distance,2D NN cluster distance,ID
            for j=1:nNNs
                remainId=remainIds(groupNNIdx(i,j));
                line1Ends=[groupFeatures(toMergeId,4),groupFeatures(toMergeId,4)+groupFeatures(toMergeId,5)];
                line2Ends=[groupFeatures(remainId,4),groupFeatures(remainId,4)+groupFeatures(remainId,5)];
                lineSegs=[(line2Ends(2)-line1Ends(1)),(line1Ends(2)-line2Ends(1))];
                verticalOverlapRatio=min(lineSegs)/max(lineSegs);
                if numel(groupHulls{toMergeId})>3 && numel(groupHulls{remainId})>3
                    horizontalOverlapRatio=overlapping(groupHulls{toMergeId},groupHulls{remainId});
                else
                    horizontalOverlapRatio=0.0;
                end
                nnClusterCentroids=clusterGroupMap(vertcat(clusterMapVGroup{remainId}),:);

                [idx,nnDs]=jitknnsearch(nnClusterCentroids(:,1:3),currentClusterCentroids(:,1:3));
                min3DSpacing=min(nnDs);

                min2DSpacing=norm(mean(nnClusterCentroids(:,1:2),1)-mean(currentClusterCentroids(:,1:2),1));
                filterMetrics(j,:)=[horizontalOverlapRatio,verticalOverlapRatio,min3DSpacing,min2DSpacing,remainId];
            end
            filterMetrics(filterMetrics(:,2)<=0,2)=0;

%             score=(filter_metrics(:,1)+PR_DampenVerticalWeight*filter_metrics(:,2))./(min(filter_metrics(:,3),filter_metrics(:,4)));
            score=exp(-((1-filterMetrics(:,1))).^2-PR_VERTICAL_WEIGHT*((1-filterMetrics(:,2))).^2-((min(filterMetrics(:,3),filterMetrics(:,4))/sigmaD).^2));
 

            [scoreSort,scoreSortI]=sort(score,'descend');
            if scoreSort(1)==0
               continue 
            end
            scoreSortRatio=scoreSort(1:end)/scoreSort(1);
            scoreSortCandidateIdx=find(scoreSortRatio>0.7);
            nScoreSortCandidateIdx=numel(scoreSortCandidateIdx);
            if nScoreSortCandidateIdx==1
                mergeNNId=groupU(filterMetrics(scoreSortI(scoreSortCandidateIdx(1)),end));
            elseif nScoreSortCandidateIdx > 1
                [filterMinSpacingD,filterMinSpacingIdx]=min(filterMetrics(scoreSortI(scoreSortCandidateIdx),3));
                mergeNNId=groupU(filterMetrics(scoreSortI(filterMinSpacingIdx),end));
            else
                continue
            end
            segs(groupVGroup{toMergeIds(i)},end)=mergeNNId;
        end
    
        clusterGroupMap(:,end)=segs(clusterUIdx,end);
        [clusterMapU,~,clusterMapV] = unique(clusterGroupMap(:,end));
        [clusterMapVSorted,clusterMapVIdx]=sort(clusterMapV);
        clusterMapVGroup = mat2cell(clusterMapVIdx, diff([0;find([(diff(clusterMapVSorted) ~= 0);1])]),1);
        
        iter=iter+1;
        
        mergedRemainIds=[mergedRemainIds;groupU(remainIds,end)];
    end
    
    unmergeIds=setdiff(groupU,mergedRemainIds);
    mergedRemainIds=unique(mergedRemainIds);
    
    unmergeIds=find(ismember(clusterMapU,unmergeIds));
    mergedRemainIds=find(ismember(clusterMapU,mergedRemainIds));
    
    [groupNNIdx,~] = jitknnsearch(groupFeatures(unmergeIds,1:2),groupFeatures(mergedRemainIds,1:2),min(PR_MIN_NN3,numel(mergedRemainIds)));
    
    nUnmergeIds=numel(unmergeIds);
    nNNs=size(groupNNIdx,2);
    for i=1:nUnmergeIds
        unmergeId=unmergeIds(i);
        currentClusterCentroids=clusterGroupMap(vertcat(clusterMapVGroup{unmergeId}),:);        
        filterMetrics=zeros(nNNs,2);%3D NN cluster distance
        for j=1:nNNs
            mergedRemainId=mergedRemainIds(groupNNIdx(i,j));
            nnClusterCentroids=clusterGroupMap(vertcat(clusterMapVGroup{mergedRemainId}),:);

            [idx,nnDs]=knnsearch(currentClusterCentroids(:,1:3),nnClusterCentroids(:,1:3));
            min3DSpacing=min(nnDs);
            filterMetrics(j,:)=[min3DSpacing,mergedRemainId];
        end
        [filterMinSpacingD,filterMinSpacingIdx]=min(filterMetrics(:,1));
        mergeNNId=groupU(filterMetrics(filterMinSpacingIdx,end));
        segs(groupVGroup{unmergeIds(i)},end)=mergeNNId;
    end

    if sum(segs(:,end))==0
        segs(:,end)=segs(:,end-1);
    end
    
    %% 4. Export laz file
    predFile=strcat(foutpath,'\',fname,'_segs.laz');
    newHdr(1) = struct(...
    'name', 'init_segs', ...
    'type', 6, ... 
    'description', 'init segs', ...
    'scale', 1.0, ...
    'offset', 0 ...
    );
    newHdr(2) = struct(...
    'name', 'interim_segs', ...
    'type', 6, ... 
    'description', 'interim segs', ...
    'scale', 1.0, ...
    'offset', 0 ...
    );
    newHdr(3) = struct(...
    'name', 'segs', ...
    'type', 6, ... 
    'description', 'segs', ...
    'scale', 1.0, ...
    'offset', 0 ...
    );

    sOut=s;
    if isfield(hdr,'attributes')
        exHdr=[hdr.attributes,newHdr];
        sOut.attribute_info = exHdr;
        sOut.attributes(:,end+1:end+3) = segs(:,end-2:end);
    else
        sOut.attribute_info = newHdr;
        sOut.attributes=[];
        sOut.attributes(:,end+1:end+3) = segs(:,end-2:end);
    end


%     mat2las(s_out, ['-o ',pred_file]);
    mat2las(sOut,[' -o ',predFile, sprintf(' -fgi_scale %.9f %.9f %.9f',hdr.x_scale_factor,hdr.y_scale_factor,hdr.z_scale_factor),  sprintf(' -fgi_offset %.9f %.9f %.9f',hdr.x_offset,hdr.y_offset,hdr.z_offset)]);
    ii=1;

    fprintf(' @%s Done\n',datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
end
