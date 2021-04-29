clear all
close all
clc

% Loading saved data
fileName = "InfMat9.mat";
load(fileName)

% Creates masks of electrodes
MirrorMask = zeros(MirrorGridSize,MirrorGridSize);
ElectMasks = zeros(MirrorGridSize,MirrorGridSize,ElectGrid^2);
for j=1:MirrorGridSize
    for i=1:MirrorGridSize
        for k = 1:size(ElectCorners,1)
            if(    ElectCorners(k,1)<=mirrorXGrid(i,j) ...
                    && mirrorXGrid(i,j)<=ElectCorners(k,2) ...
                && ElectCorners(k,3)<=mirrorYGrid(i,j) ...
                    && mirrorYGrid(i,j)<=ElectCorners(k,4) )
                ElectMasks(i,j,k) = 1;
            end
        end
        if( ~isnan(InfFuncs(i,j,1)) )
            MirrorMask(i,j) = 1;
        end
    end
end

% Accommodates the values of the influence functions to be multiplied by
% the vectors such that the defelction is z = MCal * P with P being the
% vector of preassure of each electrode
totalZPoints = sum(sum(MirrorMask));
AvgElectMask = NaN(ElectGrid^2,totalZPoints);
MCal = NaN(totalZPoints,ElectGrid^2);
for k=1:ElectGrid^2
    AvgElectMask(k,:) = ( MatUtils.matrixToVecIdxMap(ElectMasks(:,:,k),MirrorMask)...
                                    / sum(sum( ElectMasks(:,:,k) )) )';
    [vec,mirrorMaskIdxMap] = MatUtils.matrixToVecIdxMap(InfFuncs(:,:,k),MirrorMask);
    MCal(:,k) = vec;
end

K = diag(ones(ElectGrid^2,1));
MMCal = eye(totalZPoints,totalZPoints) + MCal*K*AvgElectMask;
MirrorMat = inv(MMCal)*MCal;

% YPhi = ones(ElectGrid^2,1)*0.5;
YPhi = rand(ElectGrid^2,1)*0.5 - 0.25
% Z = inv(MMCal'*MMCal)*MMCal'*MCal*YPhi;
Z = MirrorMat*YPhi;

ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
surf(mirrorXGrid,mirrorYGrid,ZMat); hold off;