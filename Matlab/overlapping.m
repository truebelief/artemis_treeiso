function overlapRatio=overlapping(convHull1,convHull2)
%OVERLAPPING Calculate the 2D horizontal overlapping ratio between two convex hulls
%   overlap_ratio = OVERLAPPING_PAIR(X,Y) 
%   Project convex hull points into binary images and calculate the maximum
%   overlapping ratio between the intersection and either X or Y
%
%   Example 1:
%   ---------
%       X = [7.23,7.46;7.32,8.31;7.31,8.34;6.68,8.67;6.64,8.69;6.62,8.68;6,8.21;5.86,8.03;5.85,8.01;5.47,6.93;5.46,6.69;5.70,6.18;5.86,6.02;6,6.04;6.11,6.07;6.13,6.08;6.60,6.39;6.63,6.41;6.90,6.68;7.23,7.46];
%       Y = [2.93,7.94;3.50,5.56;3.80,4.94;3.82,4.92;3.85,4.90;6.53,4.74;6.59,4.75;6.80,4.82;7.12,5.66;7.21,5.90;7.21,5.91;6.66,8.31;6.53,8.58;5.89,9.76;5.74,9.96;5.48,9.99;5.12,9.94;4.16,9.60;3.83,9.45;3.76,9.41;3.73,9.38;3.27,8.91;2.93,7.94];
%       ratio = overlapping(X,Y);
%   Copyright 2022 Zhouxin Xi (zhouxin.xi@uleth.ca).

scale=10;%arbitary, just used to enlarge the mask image
convHullCombo=[convHull1;convHull2];
digitizeSize = double(ceil(max(convHullCombo(:,1:2),[],1)*scale)+10);
digitizeSize=[digitizeSize(2),digitizeSize(1)];

convMask1 = poly2mask(convHull1(:,1)*scale, convHull1(:,2)*scale, digitizeSize(1),digitizeSize(2));
convMask2 = poly2mask(convHull2(:,1)*scale, convHull2(:,2)*scale, digitizeSize(1),digitizeSize(2));
% conv_union_mask= conv_mask1 | conv_mask2;
convIntersectMask = convMask1 & convMask2;
% figure,imagesc(conv_intersect_mask);
% axis equal;

% overlap_ratio=sum(conv_intersect_mask>0)/sum(conv_union_mask>0);
overlapRatio=max(sum(convIntersectMask>0)/sum(convMask1>0),sum(convIntersectMask>0)/sum(convMask2>0));
end