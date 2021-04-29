close all
clc

% Deflection axis limits
zMax = 3.5*max(max(MirrorMat));
% zMax = 0.7;
zMin = -0.2;

% Creates and set config for axes handler and surf handler
videoName = 'SimpleVideo';
vidWriter = VideoWriter(videoName,'MPEG-4');
isRecording = false;
axHand = axes;
surfHand = surf( mirrorXGrid,...
                 mirrorYGrid,...
                 NaN(MirrorGridSize,MirrorGridSize),...
                 'edgecolor','interp' );
view([-32,40])
% [axHand.View(1),axHand.View(2)]
caxis([zMin,zMax]);
colorbar; colormap jet;
shading interp
daspect([1,1,1.33]); % axis equal
zlim([zMin,zMax])

xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');
drawnow

if(isRecording)
    open(vidWriter);
    frame = getframe(gcf);
    writeVideo(vidWriter,frame);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% steps = 50;
% YPhi = zeros(ElectGrid^2,1);
% YPhiInterp = [linspace(0,1,steps)';linspace(1,0,steps)'];
% for k=1:ElectGrid^2
%     for i=1:length(YPhiInterp)
%         YPhi(k,1) = YPhiInterp(i);
%         Z = MirrorMat*YPhi;
%         ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,...
%                                 MirrorGridSize,MirrorGridSize,NaN);
%         set(surfHand,'ZData',ZMat);
%         zlim([zMin,zMax]);
%         drawnow
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steps = 80;

applyInput(PhiArr,inputMin,1,ones(ElectGrid^2,ElectGrid^2));
applyInput(PhiArr,inputMax,1,2*rand(ElectGrid^2,ElectGrid^2));
applyInput(PhiArr,0,1,ones(ElectGrid^2,ElectGrid^2));

voltageInput = [linspace(0,inputMax,steps)';...
    linspace(inputMax,-800,steps)';...
    linspace(-800,0,steps)'];
for i=1:length(voltageInput)
    applyInput(PhiArr,voltageInput(i),1,ones(ElectGrid^2,ElectGrid^2));
    YPhi = getPhiArrOutputs(PhiArr)/1000;

    Z = MirrorMat*YPhi;
    ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,...
                            MirrorGridSize,MirrorGridSize,NaN);
    set(surfHand,'ZData',ZMat);
    zlim([zMin,zMax]);
%     drawnow limitrate
    
    if(isRecording)
        drawnow
        frame = getframe(gcf);
        writeVideo(vidWriter,frame);
    end
end

steps = 40;
voltageInput = [linspace(0,inputMax*1,steps)';...
    linspace(inputMax*1,0,steps)'];
for k=1:ElectGrid^2
    for i=1:length(voltageInput)
        applyInput(PhiArr,voltageInput(i),k,ECoup);
        YPhi = getPhiArrOutputs(PhiArr)/1000;
        
        Z = MirrorMat*YPhi;
        ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,...
                                MirrorGridSize,MirrorGridSize,NaN);
        set(surfHand,'ZData',ZMat);
        zlim([zMin,zMax]);
%         drawnow limitrate
        
        if(isRecording)
            drawnow
            frame = getframe(gcf);
            writeVideo(vidWriter,frame);
        end
    end
end

if(isRecording)
    close(vidWriter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function applyInput(PhiArr,input,muxCh,ECoup)
    inputs = ECoup(muxCh,:)*input;
    for i=1:length(PhiArr)
        PhiArr(i).updateRelays(inputs(i));
    end
end

function YPhi = getPhiArrOutputs(PhiArr)
    YPhi = zeros(length(PhiArr),1);
    for i=1:length(PhiArr)
        YPhi(i) = PhiArr(i).getOutput();
    end
end

