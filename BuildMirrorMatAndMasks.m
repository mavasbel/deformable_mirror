% clear all
close all
clc

% Loading saved data
fileName = "InfMat25BigElect.mat";
load(fileName)

% ZernikePol
nZer = 4; mZer = 2;
ZenikeEfectSurfFact = 1.0;
ZernikeAmpFact = 0.02;
ZernikeRadReductFact = 1;
ZernikePol = NaN(MirrorGridSize,MirrorGridSize);
ZernikePol(MirrorMask)=0;
ZernikeMask = NaN(MirrorGridSize,MirrorGridSize);

% Creates masks of electrodes
% MirrorMask = zeros(MirrorGridSize,MirrorGridSize);
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
        r = (sqrt(mirrorXGrid(i,j)^2 + mirrorYGrid(i,j)^2)...
                /(EfectiveSurfFact*MaxRad))...
                /ZenikeEfectSurfFact;
        if( MirrorMask(i,j) && r<ZernikeRadReductFact)
            phi = atan2(mirrorYGrid(i,j),mirrorXGrid(i,j));
            ZernikePol(i,j) = zernfun(nZer,mZer,r,phi,'norm');
            ZernikeMask(i,j) = 1;
        else
            ZernikeMask(i,j) = 0;
        end
        
    end
end

% Plot Zernike function
ZernikePolMax = max(max(ZernikePol)); ZernikePolMin = min(min(ZernikePol));
ZernikePolAmp = ZernikePolMax - ZernikePolMin;
ZernikePol = 2*(ZernikePol/ZernikePolAmp)*ZernikeAmpFact;

ZernikePolMax = max(max(ZernikePol)); ZernikePolMin = min(min(ZernikePol));
ZernikePolAmp = ZernikePolMax - ZernikePolMin;
ZernikePol = ZernikePol -ZernikePolMin;

ZernikePolMax = max(max(ZernikePol)); ZernikePolMin = min(min(ZernikePol));
ZernikePolAmp = ZernikePolMax - ZernikePolMin;

figure; surf(mirrorXGrid,mirrorYGrid,ZernikePol); hold off; drawnow;
title("Zernike: n="+nZer+", m="+mZer);
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
zlim([ZernikePolMin-ZernikePolAmp*0.1,ZernikePolMax+ZernikePolAmp*0.1]);
[ZernikePolVec,~] = MatUtils.matrixToVecIdxMap(ZernikePol,MirrorMask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything here is computed for whole mirror surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accommodates the values of the influence functions to be multiplied by
% the vectors such that the defelction is z = MCal * P with P being the
% vector of preassure of each electrode
maskToUse = MirrorMask;
TotalZPoints = sum(sum(maskToUse));
AvgElectVecMask = NaN(TotalZPoints,ElectGrid^2);
ElectVecMask = zeros(TotalZPoints,ElectGrid);
MCal = NaN(TotalZPoints,ElectGrid^2);
for k=1:ElectGrid^2
    [ElectVecMask(:,k),~] = MatUtils.matrixToVecIdxMap(ElectMasks(:,:,k),maskToUse);
    AvgElectVecMask(:,k) = ElectVecMask(:,k)/sum(ElectVecMask(:,k));
    [vec,mirrorMaskIdxMap] = MatUtils.matrixToVecIdxMap(InfFuncs(:,:,k),maskToUse);
    MCal(:,k) = vec;
end

% Computes spring constans matrix K and final mirror mat such that
% Z = MirrorMat * YPhi with
% MirrorMat = inv(I+Matcal*K*ElectAvg)*MCal such t
K = diag(ones(ElectGrid^2,1));
MMCal = eye(TotalZPoints,TotalZPoints) + MCal*K*(AvgElectVecMask');
MirrorMat = inv(MMCal)*MCal;
pinvMirrorMat = inv(MirrorMat'*MirrorMat)*MirrorMat';

% Plots illustrative output with some random pressures
% YPhi = ones(ElectGrid^2,1)*0.5;
YPhi = rand(ElectGrid^2,1)*0.5 - 0.25;
% Z = inv(MMCal'*MMCal)*MMCal'*MCal*YPhi;
Z = MirrorMat*YPhi;
ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
figure; surf(mirrorXGrid,mirrorYGrid,ZMat); hold off; drawnow;
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
zlim([-1,1]*1.25);
title("Random Pressures");

% Some Influence functions
figure; counter=1;
for i=[1 5 10 15 20 25] %i=1:ElectGrid^2
    handler = subplot(2,3,counter);
    YPhi = [zeros(i-1,1);1;zeros(ElectGrid^2-i,1)];
    Z = MirrorMat*YPhi;
    ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
    surf(handler,mirrorXGrid,mirrorYGrid,ZMat,'edgecolor','interp'); hold off; drawnow;
    xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
    zlim([-1,1]*0.30);
    title("Actuator "+i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything here is computed using zernike mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accommodates the values of the influence functions to be multiplied by
% the vectors such that the defelction is z = MCal * P with P being the
% vector of preassure of each electrode
maskToUse = ZernikeMask;
ZerTotalZPoints = sum(sum(maskToUse));
ZerAvgElectVecMask = NaN(ZerTotalZPoints,ElectGrid^2);
ZerElectVecMask = zeros(ZerTotalZPoints,ElectGrid);
ZerMCal = NaN(ZerTotalZPoints,ElectGrid^2);
for k=1:ElectGrid^2
    [ZerElectVecMask(:,k),~] = MatUtils.matrixToVecIdxMap(ElectMasks(:,:,k),maskToUse);
    ZerAvgElectVecMask(:,k) = ZerElectVecMask(:,k)/sum(ZerElectVecMask(:,k));
    [vec,~] = MatUtils.matrixToVecIdxMap(InfFuncs(:,:,k),maskToUse);
    ZerMCal(:,k) = vec;
end

% Computes spring constans matrix K and final mirror mat such that
% Z = MirrorMat * YPhi with
% MirrorMat = inv(I+Matcal*K*ElectAvg)*MCal such t
K = diag(ones(ElectGrid^2,1));
ZerMMCal = eye(ZerTotalZPoints,ZerTotalZPoints) + ZerMCal*K*(ZerAvgElectVecMask');
ZerMirrorMat = inv(ZerMMCal)*ZerMCal;
ZerPinvMirrorMat = inv(ZerMirrorMat'*ZerMirrorMat)*ZerMirrorMat';

% Compute reference values with MLSE
[ZernikePolMaskVec,ZernikePolMaskVecIdxMap] = MatUtils.matrixToVecIdxMap(ZernikePol,maskToUse);
YPhiRef = ZerPinvMirrorMat*ZernikePolMaskVec;
ZFit = MirrorMat*YPhiRef;
ZFitMat = MatUtils.vecIdxMapToMatrix(ZFit,mirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
ZFitZerMask = MatUtils.matrixToVecIdxMap(ZFitMat,maskToUse);

figure; surf(mirrorXGrid,mirrorYGrid,ZFitMat,'edgecolor','interp'); hold off; drawnow;
colorbar; colormap jet;
shading interp
ZMax = max(ZFit); ZMin = min(ZFit); ZPad = 0.20;
ZAmp = ZMax - ZMin; ZMaxAbs = max(abs([ZMin,ZMax]));
zlim([-1 1]*ZMaxAbs); caxis([-1,1]*ZMaxAbs);
daspect([1,1,max([ZAmp,ZPad])*1.35]);
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]); 
title("Reference $Z_d$",'interpreter','latex');
% zlim([ZernikePolMin-ZernikePolAmp*0.1, ZernikePolMax+ZernikePolAmp*0.1]);
% daspect([1,1,0.05]);
% title("Best Zernike Fit");
TotalSqrError = sum((ZFit-ZernikePolVec).^2)

