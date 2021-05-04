close all
clc

% Deflection axis limits
ZMax = 3.4*max(max(MirrorMat)); ZMin = -0.75*max(max(MirrorMat));
ZAmp = ZMax - ZMin;
ZPad = 0.20;

% Creates and set config for axes handler and surf handler
if(exist('vidWriter')) close(vidWriter); end
videoName = 'ControlTestVideo';
vidWriter = VideoWriter(videoName,'MPEG-4');
isRecording = false;
axHand = axes;
surfHand = surf( mirrorXGrid,...
                 mirrorYGrid,...
                 NaN(MirrorGridSize,MirrorGridSize),...
                 'edgecolor','interp' );
view([-32,40])
% [axHand.View(1),axHand.View(2)]
colorbar; colormap jet;
shading interp
caxis([ZMin-ZAmp*ZPad, ZMax+ZAmp*ZPad]);
zlim([ZMin-ZAmp*ZPad, ZMax+ZAmp*ZPad]);
daspect([1,1,1.33]); % axis equal

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
    close(vidWriter)
    open(vidWriter);
    frame = getframe(gcf);
    writeVideo(vidWriter,frame);
end

% Parameters for plot and control
ZFitMax = max(ZFit); ZFitMin = min(ZFit);
AmpZFit = ZFitMax - ZFitMin;
PadZFit = 0.1;
samplesPerSegment = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the mirror resetting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title("Resetting");
voltageInput = linspace(0,inputMax,samplesPerSegment);
[Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(inputMax,0,samplesPerSegment);
[Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(0,-800,samplesPerSegment);
[Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(-800,0,samplesPerSegment);
[Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1200;
iter = 0;

[ZZerMask,~] = MatUtils.matrixToVecIdxMap(ZMat,ZernikeMask);
ZError = ZFitZerMask - ZZerMask;
PhiError = ZerPinvMirrorMat*ZError;
ElectInputs = ones(ElectGrid^2,1)*600;

TotalZSquareError = sum(ZError.^2);

while(true)
    disp("-------------------------------");
    disp("Iteration: " + (iter+1));
    disp("Deflection Total Square Error: " + sum(ZError.^2));
    disp("Pressure Total Square Error: " + sum(PhiError.^2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     ZError = ZFit - Z;
%     PhiError = pinvMirrorMat*ZError;    
    [ZZerMask,~] = MatUtils.matrixToVecIdxMap(ZMat,ZernikeMask);
    ZError = ZFitZerMask - ZZerMask; 
    TotalZSquareError = sum(ZError.^2);
    PhiError = ZerPinvMirrorMat*ZError;
    
    if TotalZSquareError<=(AmpZFit/100)*5 break; end

    ElectInputs = ElectInputs + lambda*PhiError;
    [sortedInputs,idx] = sort(ElectInputs,'descend');    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ResetFlag = false;
    for i=1:length(PhiError)
        if PhiError(i)<0 ResetFlag = true; break; end
    end
    
    if(ResetFlag)
        disp("Resetting...");
        title("Resetting");
        voltageInput = linspace(0,inputMax,samplesPerSegment);
        [Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
        voltageInput = linspace(inputMax,0,samplesPerSegment);
        [Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
        voltageInput = linspace(0,-800,samplesPerSegment);
        [Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
        voltageInput = linspace(-800,0,samplesPerSegment);
        [Z,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:length(sortedInputs)
        title("Iteration: " + (iter+1) + ", Actuating Electrode: " + idx(i));
        disp("Actuating step: " + i);
        voltageInput = [linspace(0,sortedInputs(i),samplesPerSegment)';...
            linspace(sortedInputs(i),0,samplesPerSegment)'];
        [Z,ZMat,Phi] = applyInputAnim(voltageInput,idx(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iter = iter + 1;
end

disp("-------------------------------");
disp("Total Iteration: " + iter);
disp("Final Deflection Total Square Error: " + sum(ZError.^2));
disp("Fintal Pressure Total Square Error: " + sum(PhiError.^2));

if(isRecording)
    close(vidWriter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z,ZMat,YPhi] = applyInputAnim(input,elect,forceReload)
    persistent PhiArr ECoup MirrorMat mirrorMaskIdxMap MirrorGridSize ...
                surfHand isRecording vidWriter ZPad;
    if isempty(PhiArr) || (nargin==3 && forceReload)
        PhiArr = evalin('base','PhiArr');
        ECoup = evalin('base','ECoup');
        MirrorMat = evalin('base','MirrorMat');
        mirrorMaskIdxMap = evalin('base','mirrorMaskIdxMap');
        MirrorGridSize = evalin('base','MirrorGridSize');
        surfHand = evalin('base','surfHand');
        isRecording = evalin('base','isRecording');
        vidWriter = evalin('base','vidWriter');
        ZPad = evalin('base','ZPad');
    end
    
    for i=1:length(input)
        % Apply input
        applyInput(PhiArr,input(i),1,ECoup(elect,:));
        YPhi = getPhiArrOutputs(PhiArr);

        % Compute deflections
        Z = MirrorMat*YPhi;
        ZMat = MatUtils.vecIdxMapToMatrix(Z,mirrorMaskIdxMap,...
                                MirrorGridSize,MirrorGridSize,NaN);

        % Config limits and color
        ZMax = max(Z); ZMin = min(Z); 
        ZAmp = ZMax - ZMin; ZMaxAbs = max(abs([ZMin,ZMax]));
        zlim([-1 1]*ZMaxAbs); caxis([-1,1]*ZMaxAbs);
        daspect([1,1,max([ZAmp,ZPad])*1.35]);
        
        % Update plot data
        set(surfHand,'ZData',ZMat);         
        
        % Record
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

function Phi = getPhiArrOutputs(PhiArr)
    Phi = zeros(length(PhiArr),1);
    for i=1:length(PhiArr)
        Phi(i) = PhiArr(i).getOutput();
    end
    Phi = Phi/1000; % This is a scale factor
end
