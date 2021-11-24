% clear all
close all
clc

% Loading saved data
% fileName = "InfMat25BigElect.mat";
fileName = "InfMat_5x5_Eff_55_Wid_12.mat";
load(fileName)

% Influence functions scaled by the relation between radius and surface tension
% ScaledInfFuncs = InfFuncs/max(max(max(InfFuncs)))-min(min(min(InfFuncs))) ); % Assuming surface tension and radious relation is rad = 10cm, tension=100N
ScaledInfFuncs = InfFuncs; % Assuming surface tension and radious relation is rad = 10cm = 100mm, tension=10N

% ZernikePol
ZerN = 2; ZerM = -2;
% ZerN = 1; ZerM = -1;
% ZerN = 4; ZerM = 2;
ZenEffectSurfFact = 1/0.60;
% ZerAmpFact = 0.025;
ZerAmpFact = 50;
ZerRadTrim = 0.65;
ZerPol = NaN(MirrorGridSize,MirrorGridSize);
ZerMask = zeros(MirrorGridSize,MirrorGridSize);

% Creates masks of electrodes
% MirrorMask = zeros(MirrorGridSize,MirrorGridSize);
ElectMasks = zeros(MirrorGridSize,MirrorGridSize,ElectGrid^2);

% Compute Zernike polynomial and its mask
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
                /ZenEffectSurfFact;
        if( MirrorMask(i,j) && r<ZerRadTrim)
            phi = atan2(mirrorYGrid(i,j),mirrorXGrid(i,j));
            ZerPol(i,j) = zernfun(ZerN,ZerM,r,phi,'norm');
            ZerMask(i,j) = 1;
        end
        
    end
end

% % Multiplot Zerenikes
% ZerPolNs = [1 1 2 2 2 3]; ZerPolMs = [-1 1 -2 0 2 1];
% ZerPols = NaN(MirrorGridSize,MirrorGridSize,length(ZerPolNs));
% figure; counter=1;
% for k = 1:length(ZerPolNs)
%     for j=1:MirrorGridSize
%         for i=1:MirrorGridSize
%             r = sqrt(mirrorXGrid(i,j)^2 + mirrorYGrid(i,j)^2);
%             if( r<0.95 )
%                 phi = atan2(mirrorYGrid(i,j),mirrorXGrid(i,j));
%                 ZerPols(i,j,k) = zernfun(ZerPolNs(k),ZerPolMs(k),r,phi,'norm');
%             end
%         end
%     end
%     handler = subplot(2,3,counter);
%     surf(handler,mirrorXGrid,mirrorYGrid,ZerPols(:,:,k),'edgecolor','interp'); hold off;
%     colormap jet; shading interp;
%     xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
%     zlim([-1,1]*1.5);
%     title("Zernike n: " + ZerPolNs(k) + ", m: " + ZerPolMs(k)); drawnow;
%     counter = counter + 1;
% end

% Plot Zernike polynomial
ZerPolMax = max(max(ZerPol)); ZerPolMin = min(min(ZerPol));
ZerPolAmp = ZerPolMax - ZerPolMin;
ZerPol = 2*(ZerPol/ZerPolAmp)*ZerAmpFact;

ZerPolMax = max(max(ZerPol)); ZerPolMin = min(min(ZerPol));
ZerPolAmp = ZerPolMax - ZerPolMin;
ZerPol = ZerPol - ZerPolMin;

ZerPolMax = max(max(ZerPol)); ZerPolMin = min(min(ZerPol));
ZerPolAmp = ZerPolMax - ZerPolMin;

figure; surf(mirrorXGrid,mirrorYGrid,ZerPol); hold off;
colormap jet; shading interp; colorbar;
title("Zernike: n="+ZerN+", m="+ZerM);
view([-32,40])
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
% zlim([ZerPolMin-ZerPolAmp*0.1,ZerPolMax+ZerPolAmp*0.1]);
zlim([-80,180]); caxis([-30,120]);
daspect([1,1,3*120]); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything here is computed for whole mirror surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accommodates the values of the influence functions to be multiplied by
% the vectors such that the defelction is z = MCal * P with P being the
% vector with the preassure of each electrode
maskToUse = MirrorMask;
ZPointsTotal = sum(sum(maskToUse));
AvgElectMaskVec = NaN(ZPointsTotal,ElectGrid^2);
ElectMaskVec = zeros(ZPointsTotal,ElectGrid);
ElectAreaPercent = (sum(sum(sum(ElectMasks)))/ElectGrid^2)/ZPointsTotal;
MCal = NaN(ZPointsTotal,ElectGrid^2);
for k=1:ElectGrid^2
    [ElectMaskVec(:,k),MirrorMaskIdxMap] = MatUtils.matrixToVecIdxMap(ElectMasks(:,:,k),maskToUse);
    AvgElectMaskVec(:,k) = ElectMaskVec(:,k)/sum(ElectMaskVec(:,k));
    [vec,~] = MatUtils.matrixToVecIdxMap(ScaledInfFuncs(:,:,k),maskToUse);
    MCal(:,k) = vec;
end

% Creates spring constans matrix K and computes HDM matrix H
K = 170*diag(ones(ElectGrid^2,1))/(sum(sum(sum(ElectMasks)))/ElectGrid^2) ; % Stiffness 170 N/m
HBold = inv( eye(ZPointsTotal,ZPointsTotal)+MCal*K*(AvgElectMaskVec') )*MCal;
pInvHBold = inv(HBold'*HBold)*HBold';

% % Some Influence functions
% figure; counter=1;
% for i=[1 5 10 15 20 25] %i=1:ElectGrid^2
%     handler = subplot(2,3,counter);
%     ZMat = MatUtils.vecIdxMapToMatrix(HBold(:,i),MirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
%     surf(handler,mirrorXGrid,mirrorYGrid,ZMat,'edgecolor','interp'); hold off;
%     colormap jet; shading interp;
%     xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
%     zlim([ZerPolMin-ZerPolAmp*0.1,ZerPolMax+ZerPolAmp*0.1]);
%     title("Actuator "+i); drawnow;
%     counter = counter + 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything here is computed using Zernike mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maskToUse = ZerMask;
ZerZPointsTotal = sum(sum(maskToUse));
ZerAvgElectMaskVec = NaN(ZerZPointsTotal,ElectGrid^2);
ZerElectMaskVec = zeros(ZerZPointsTotal,ElectGrid);
ZerMCal = NaN(ZerZPointsTotal,ElectGrid^2);
for k=1:ElectGrid^2
    [ZerElectMaskVec(:,k),ZerMaskIdxMap] = MatUtils.matrixToVecIdxMap(ElectMasks(:,:,k),maskToUse);
    ZerAvgElectMaskVec(:,k) = ZerElectMaskVec(:,k)/sum(ElectMaskVec(:,k));
    [vec,~] = MatUtils.matrixToVecIdxMap(ScaledInfFuncs(:,:,k),maskToUse);
    ZerMCal(:,k) = vec;
end

% Creates spring constans matrix K and computes HDM matrix H
ZerHBold = inv( eye(ZerZPointsTotal,ZerZPointsTotal)+ZerMCal*K*(ZerAvgElectMaskVec') )*ZerMCal;
ZerPInvHBold = inv(ZerHBold'*ZerHBold)*ZerHBold';

% Estimate best pressure values with MLSE for testing
[ZerPolVec_ZerMask,~] = MatUtils.matrixToVecIdxMap(ZerPol,maskToUse);
PressRef = ZerPInvHBold*ZerPolVec_ZerMask;
ZFit = HBold*PressRef;
% sum(abs(ZerPInvHBold*ZerPolVec_ZerMask - pInvHBold*ZFit)) % This compares Press[using Zernike mask] and Press[using mirror mask]
ZFitMat = MatUtils.vecIdxMapToMatrix(ZFit,MirrorMaskIdxMap,MirrorGridSize,MirrorGridSize,NaN);
ZFitVecZerMask = MatUtils.matrixToVecIdxMap(ZFitMat,maskToUse);
ZFitMax = max(ZFit); ZFitMin = min(ZFit); ZFitPad = 0.20;
ZFitAmp = ZFitMax - ZFitMin; ZFitMaxAbs = 100;%max(abs([ZFitMin,ZFitMax]));

% Plot error between Zernike polynomial and best approximation
figure; surf(mirrorXGrid,mirrorYGrid,(ZerPol-ZFitMat).*ZerMask,'edgecolor','interp'); hold off;
colormap jet; shading interp; colorbar;
title('ZerPol - ZFit');
view([-32,40])
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]); 
zlim([-80,180]); caxis([-30,120]);
daspect([1,1,3*120]); drawnow;

% Error
% TotalError = nansum(nansum( ((ZerPol-ZFitMat).*ZerMask).^2 ))/ZPointsTotal
MaxError = max(max(ZerPol-ZFitMat))

% Plot best approximation
figure; surf(mirrorXGrid,mirrorYGrid,ZFitMat,'edgecolor','interp'); hold off;
colormap jet; shading interp; colorbar;
% title('ZFit');
view([-32,40])
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]);
zlim([-80,180]); caxis([-30,120]);
daspect([1,1,3*120]); drawnow;

% Labels for saving the figure
xlabel('','interpreter','latex');
ylabel('','interpreter','latex');
zlabel('$nm$','interpreter','latex'); drawnow;
saveas(gca, "hdm_mirror_reference","epsc");

