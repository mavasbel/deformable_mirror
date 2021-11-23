close all
clc

% Deflection axis limits
zMax = 3.5*max(max(HBold));
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
colorbar; colormap jet; shading interp;
daspect([1,1,1.33]); % axis equal
zlim([zMin,zMax])

% Set axis labels
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');

% Initialize mirror with random remnants
applyInput(PhiArr,inputMin,1,ones(ElectGrid^2,ElectGrid^2));
applyInput(PhiArr,inputMax,1,2*rand(ElectGrid^2,ElectGrid^2));
applyInput(PhiArr,0,1,ones(ElectGrid^2,ElectGrid^2));
drawnow;
if(isRecording)
    open(vidWriter);
    frame = getframe(gcf);
    writeVideo(vidWriter,frame);
end

steps = 80;
voltageInput = [linspace(0,inputMax,steps)';...
    linspace(inputMax,-800,steps)';...
    linspace(-800,0,steps)'];
applyInputAnim(voltageInput,ElectGrid^2+1,true);

steps = 40;
voltageInput = [linspace(0,inputMax*1,steps)';...
    linspace(inputMax*1,0,steps)'];
for k=1:ElectGrid^2
    applyInputAnim(voltageInput,k);
end

if(isRecording)
    close(vidWriter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = applyInputAnim(input,elect,forceReload)
    persistent PhiArr ECoup HBold MirrorMaskIdxMap MirrorGridSize ...
                surfHand isRecording vidWriter;
    if isempty(PhiArr) || (nargin==3 && forceReload)
        PhiArr = evalin('base','PhiArr');
        ECoup = evalin('base','ECoup');
        HBold = evalin('base','HBold');
        MirrorMaskIdxMap = evalin('base','MirrorMaskIdxMap');
        MirrorGridSize = evalin('base','MirrorGridSize');
        surfHand = evalin('base','surfHand');
        isRecording = evalin('base','isRecording');
        vidWriter = evalin('base','vidWriter');
    end
    
    for i=1:length(input)
        applyInput(PhiArr,input(i),1,ECoup(elect,:));
        YPhi = getPhiArrOutputs(PhiArr);

        Z = HBold*YPhi;
        ZMat = MatUtils.vecIdxMapToMatrix(Z,MirrorMaskIdxMap,...
                                MirrorGridSize,MirrorGridSize,NaN);
        set(surfHand,'ZData',ZMat);

        if(isRecording)
            drawnow
            frame = getframe(gcf);
            writeVideo(vidWriter,frame);
        else
            drawnow limitrate
        end
    end
end

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
    YPhi = YPhi/1000; % This is a scale factor
end
