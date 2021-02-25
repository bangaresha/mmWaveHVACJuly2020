function [wallxyz1, wallxyz2, wallxyz3, wallxyz4,wallX,wallY,wallZ] = CSV23D_V1 (demoMode,groundLevel,ceilingLevel,varargin)
[fileName,fileAddress] = uigetfile('*.csv');
[~,~,pointData] = xlsread([fileAddress,fileName]);
for i = 1:numel(pointData)
    X(i,1) = str2num(strtok(pointData{i,1},';'));
    [~,remain{i,1}] = strtok(pointData{i,1},';');
    Y(i,1) = str2num(strtok(remain{i},';'));
end
wallxyz1 = zeros(size(X,1)/2,3);
wallxyz2 = wallxyz1;
wallxyz3 = wallxyz1;
wallxyz4 = wallxyz1;
% Defining walls clockwise
for i = 1:numel(X)/2
    wallxyz1(i,:) = [X((i*2)-1,1),Y((i*2)-1,1),groundLevel];
    wallxyz2(i,:) = [X((i*2)-1,1),Y((i*2)-1,1),ceilingLevel];
    wallxyz3(i,:) = [X(i*2,1),Y(i*2,1),ceilingLevel];
    wallxyz4(i,:) = [X(i*2,1),Y(i*2,1),groundLevel];
end
wallX = [wallxyz1(:,1)';wallxyz2(:,1)';wallxyz3(:,1)';wallxyz4(:,1)'];
wallY = [wallxyz1(:,2)';wallxyz2(:,2)';wallxyz3(:,2)';wallxyz4(:,2)'];
wallZ = [wallxyz1(:,3)';wallxyz2(:,3)';wallxyz3(:,3)';wallxyz4(:,3)'];
wallC = ones(size(wallX)); 
% figure
% fill3(wallX, wallY, wallZ,wallC)
% text(varargin{1}(1),varargin{1}(2),varargin{1}(3),'TX','Color','Red')
% title("3D")