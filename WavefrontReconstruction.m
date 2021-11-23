clear all
close all
clc

% Parameters
MaxRad = 1;
CoarseGridSize = 20;
DenseGridSize = 100;
ZerPolyMaxLevel = 7; 

%% Coarse grid computations

% Create coarse grid of points to compute Zernike polynomials
[~,~,CoarseXGrid,CoarseYGrid] = createGrid(MaxRad,CoarseGridSize);
        
% Compute first MaxZerLevel-th levels of Zernike polynomials in Coarse grid
[CoarseZerPoly,CircMask] = computeZernikes(CoarseXGrid,CoarseYGrid,ZerPolyMaxLevel);

% Stack coarse Zernike polynomials
CoarseTotalZernikes = size(CoarseZerPoly,3);
CoarseTotalZPoints = sum(sum(CircMask));
CoarseZerPolMat = zeros(CoarseTotalZPoints,CoarseTotalZernikes);
for k=1:CoarseTotalZernikes
     [vec,CoarseCircMaskIdxMap] = MatUtils.matrixToVecIdxMap(CoarseZerPoly(:,:,k), CircMask);
     CoarseZerPolMat(:,k) = vec;
end

%% Dense grid computations

% Create dense grid of points to compute Zernike polynomials
[~,~,DenseXGrid,DenseYGrid] = createGrid(MaxRad,DenseGridSize);
        
% Compute first MaxZerLevel-th levels of Zernike polynomials in Coarse grid
[DenseZerPoly,CircMask] = computeZernikes(DenseXGrid,DenseYGrid,ZerPolyMaxLevel);

% Stack dense Zernike polynomials
DenseTotalZernikes = size(DenseZerPoly,3);
DenseTotalZPoints = sum(sum(CircMask));
DenseZerPolMat = zeros(DenseTotalZPoints,DenseTotalZernikes);
for k=1:DenseTotalZernikes
     [vec,DenseCircMaskIdxMap] = MatUtils.matrixToVecIdxMap(DenseZerPoly(:,:,k), CircMask);
     DenseZerPolMat(:,k) = vec;
end

%% Random weights computations with coarse grid 

% Generate vector of random weights
randomWeights = (rand(CoarseTotalZernikes,1)-0.5)*2*5;

% Generate random vector of deformation with noise
CoarseZRandVec = CoarseZerPolMat*randomWeights + (rand(CoarseTotalZPoints,1)-0.5)*2*0.075;

% Convert to random deformations with noise to matrix and plot
CoarseZRand = MatUtils.vecIdxMapToMatrix(CoarseZRandVec,CoarseCircMaskIdxMap,CoarseGridSize,CoarseGridSize,NaN);

%% Fitting with coarse grid

identWeights = inv(CoarseZerPolMat'*CoarseZerPolMat)*CoarseZerPolMat'*CoarseZRandVec;

%% Ploting of the dense Zernike polinomials identified from coarse grid

% Generate identified vector of deformation
DenseZIdentVec = DenseZerPolMat*identWeights;

% Convert identified deformations to matrix and plot
DenseZIdent = MatUtils.vecIdxMapToMatrix(DenseZIdentVec,DenseCircMaskIdxMap,DenseGridSize,DenseGridSize,NaN);

% Plot
disp(['Error: ', num2str(sum(abs(randomWeights-identWeights).^2))])
surf(DenseXGrid,DenseYGrid,DenseZIdent,'edgecolor','interp'); hold off; drawnow;

%% Auxiliary functions

function [mirrorXIsoLine,mirrorYIsoLine,...
            mirrorXGrid,mirrorYGrid] = createGrid(MaxRad,MirrorGridSize)
    mirrorXLim = [-MaxRad, MaxRad];
    mirrorYLim = [-MaxRad, MaxRad];
    mirrorXIsoLine = linspace(mirrorXLim(1),mirrorXLim(2),MirrorGridSize);
    mirrorYIsoLine = linspace(mirrorYLim(1),mirrorYLim(2),MirrorGridSize);
    [mirrorXGrid,mirrorYGrid] = meshgrid(mirrorXIsoLine,mirrorYIsoLine);
end

function [ZernikePol,CircMask] ...
            = computeZernikes(XGrid,YGrid,MaxZerLevel)
    
    CircMask = XGrid.^2 + YGrid.^2 <= ones(size(YGrid,1),size(XGrid,2));
    ZernikePol = NaN(size(YGrid,1),size(XGrid,2));

    PolC = 1;
    for ZerN = 0:MaxZerLevel
        for ZerM = ZerN:-2:-ZerN
            for j=1:size(XGrid,2)
                parfor i=1:size(YGrid,1)
                    if( CircMask(i,j) )
                        r = (sqrt(XGrid(i,j)^2 + YGrid(i,j)^2));
                        phi = atan2(YGrid(i,j),XGrid(i,j));
                        ZernikePol(i,j,PolC) = zernfun(ZerN,ZerM,r,phi,'norm');
                    else
                        ZernikePol(i,j,PolC) = NaN;
                    end
                end
            end
            PolC = PolC+1;
        end
    end

end