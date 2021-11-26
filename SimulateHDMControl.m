close all
clc

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
             colormap jet; shading interp; colorbar;
view([-32,40])
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]); 
zlim([-80,180]); caxis([-30,120]);
daspect([1,1,3*120]); drawnow;

% Set axis labels
xlabel('','interpreter','latex');
ylabel('','interpreter','latex');
zlabel('$nm$','interpreter','latex'); drawnow;

% Initialize mirror with random remnants
[ZVec,ZMat,Press] = applyInput(inputMin,1,ones(ElectGrid^2,1),true);
[ZVec,ZMat,Press] = applyInput(inputMax,1,rand(ElectGrid^2,1));
[ZVec,ZMat,Press] = applyInput(0,1,ones(ElectGrid^2),true);
drawnow;
if(isRecording)
    close(vidWriter)
    open(vidWriter);
    frame = getframe(gcf);
    writeVideo(vidWriter,frame);
end

% Control parameters
kappa = 0.130;
ZErrorMaxThres = ZFitAmp*1/100;
samplesPerSegment = 2;
ElectInputs = 700*ones(ElectGrid^2,1); % Initial value

%Save figure
saveas(gca, "hdm_control_initial.eps","epsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the mirror resetting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title("Resetting");
voltageInput = linspace(0,inputMax,samplesPerSegment);
[ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup,true);
voltageInput = linspace(inputMax,0,samplesPerSegment);
[ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
voltageInput = linspace(0,-800,samplesPerSegment);
[ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
voltageInput = linspace(-800,0,samplesPerSegment);
[ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);

[ZVecZerMask,~] = MatUtils.matrixToVecIdxMap(ZMat,ZerMask);
ZError = ZFitVecZerMask - ZVecZerMask; 
PressError = ZerPInvHBold*ZError;
disp("Reset")
disp("Max Deflection Error Threshold: " + ZErrorMaxThres);
disp("Max Deflection Error: " + max(abs(ZError)));
disp("Max Pressure Error: " + max(abs(PressError)));
disp("Min Deflection Error: " + min(ZError));
disp("Min Pressure Error: " + min(PressError));
disp("-------------------------------");
    
%Save figure
saveas(gca, "hdm_control_reset.eps","epsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
while(true)
    
    if max(abs(ZError))<=ZErrorMaxThres break;
    else iter = iter + 1; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PressError(PressError<0) = 0;
    ElectInputs = ElectInputs + kappa*PressError;
    %[sortedInputs,idx] = sort(ElectInputs,'descend');    
    sortedInputs = ElectInputs;
    idx = [1:length(ElectInputs)]';    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     ResetFlag = false;
%     for i=1:length(PhiError)
%         if PhiError(i)<0 ResetFlag = true; break; end
%     end
%     if(ResetFlag)
%         disp("Resetting...");
%         title("Resetting");
%         voltageInput = linspace(0,inputMax,samplesPerSegment);
%         [ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
%         voltageInput = linspace(inputMax,0,samplesPerSegment);
%         [ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
%         voltageInput = linspace(0,-800,samplesPerSegment);
%         [ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
%         voltageInput = linspace(-800,0,samplesPerSegment);
%         [ZVec,ZMat,Press] = applyInputAnim(voltageInput,ElectGrid^2+1,ECoup);
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:length(sortedInputs)
%         title("Iteration: " + (iter+1) + ", Actuating Electrode: " + idx(i));
        voltageInput = [linspace(0,sortedInputs(i),samplesPerSegment)';...
            linspace(sortedInputs(i),0,samplesPerSegment)'];
        [ZVec,ZMat,Press] = applyInputAnim(voltageInput,idx(i),ECoup);
        disp("Iteration: " + iter ...
            + ", Step: " + i ...
            + ", Electrode: " + idx(i) ...
            + ", Amplitude: " + sortedInputs(i) ...
            + ", Pressure: " + Press(i));
    end
    
    [ZVecZerMask,~] = MatUtils.matrixToVecIdxMap(ZMat,ZerMask);
    ZError = ZFitVecZerMask - ZVecZerMask;
%     ZErrorPos = ZError; ZErrorPos(ZError<0)=0;
    PressError = (ZerPInvHBold*ZError);
    disp("Max Deflection Error: " + max(abs(ZError)));
    disp("Max Pressure Error: " + max(abs(PressError)));
    disp("Min Deflection Error: " + min(ZError));
    disp("Min Pressure Error: " + min(PressError));
    disp("-------------------------------");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Save figure
%     title("Iteration: " + (iter+1) + ", Actuating Electrode: " + idx(i));
    saveas(gca, "hdm_control_iter_"+iter+".eps","epsc");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

disp("Final Max Deflection Error: " + max(abs(ZError)) ...
    + " less than Max Deflection Error Threshold: " + ZErrorMaxThres);
disp("Total Iterations: " + iter);
disp("-------------------------------");

if(isRecording)
    close(vidWriter);
end

% Ploting error between best approximatio and achieved approximation
figure; surf(mirrorXGrid,mirrorYGrid,(ZFitMat-ZMat).*ZerMask,'edgecolor','interp'); hold off;
colormap jet; shading interp; colorbar;
title('ZFit - Z');
view([-32,40])
xlim([-MaxRad MaxRad]); ylim([-MaxRad MaxRad]); 
zlim([-80,180]); caxis([-30,120]);
daspect([1,1,3*120]); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ZVec,ZMat,Press] = applyInputAnim(input,muxCh,ECoup,forceReload)
    persistent surfHand isRecording vidWriter;
    if isempty(surfHand) || (nargin==4 && forceReload)
        surfHand = evalin('base','surfHand');
        isRecording = evalin('base','isRecording');
        vidWriter = evalin('base','vidWriter');
    end
    
    for i=1:length(input)
        % Apply input
        [ZVec,ZMat,Press] = applyInput(input(i),muxCh,ECoup);
        
        % Update plot data
        set(surfHand,'ZData',ZMat); drawnow;
        
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

function [ZVec,ZMat,Press] = applyInput(input,muxCh,ECoup,forceReload)
    persistent PhiArr HBold MirrorMaskIdxMap MirrorGridSize YoungMod;
    if isempty(PhiArr) || (nargin==4 && forceReload)
        PhiArr = evalin('base','PhiArr');
        HBold = evalin('base','HBold');
        MirrorMaskIdxMap = evalin('base','MirrorMaskIdxMap');
        MirrorGridSize = evalin('base','MirrorGridSize');
        YoungMod = evalin('base','YoungMod');
    end
    
    inputs = ECoup(:,muxCh)*input;
    for i=1:length(PhiArr)
        PhiArr(i).updateRelays(inputs(i));
    end
    
    % Compute deflections
    Press = YoungMod*getPhiArrOutputs(PhiArr);
    ZVec = HBold*Press;
    ZMat = MatUtils.vecIdxMapToMatrix(ZVec,MirrorMaskIdxMap,...
                            MirrorGridSize,MirrorGridSize,NaN);
end

function Phi = getPhiArrOutputs(PhiArr)
%     persistent dataHandler;
%     if isempty(dataHandler)
%         dataHandler = evalin('base','dataHandler');
%     end
    Phi = zeros(length(PhiArr),1);
    for i=1:length(PhiArr)
        Phi(i) = PhiArr(i).getOutput(); 
    end
end
