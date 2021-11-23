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
samplesPerSegment = 10;

%Save figure
saveas(gca, "hdm_control_initial.eps","epsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the mirror resetting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title("Resetting");
voltageInput = linspace(0,inputMax,samplesPerSegment);
[ZVec,ZMat,Phi] = applyCoupledInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(inputMax,0,samplesPerSegment);
[ZVec,ZMat,Phi] = applyCoupledInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(0,-800,samplesPerSegment);
[ZVec,ZMat,Phi] = applyCoupledInputAnim(voltageInput,ElectGrid^2+1,true);
voltageInput = linspace(-800,0,samplesPerSegment);
[ZVec,ZMat,Phi] = applyCoupledInputAnim(voltageInput,ElectGrid^2+1,true);

%Save figure
saveas(gca, "hdm_control_reset.eps","epsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = 0.45;
ElectInputs = ones(ElectGrid^2,1)*700;
iter = 0;
while(true)
    [ZVecZerMask,~] = MatUtils.matrixToVecIdxMap(ZMat,ZerMask);
    ZError = ZFitVecZerMask - ZVecZerMask; 
    PressError = ZerPInvHBold*ZError;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sum(ZError.^2)/MirrorGridSize^2<=ZFitAmp*0.5/100 break; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iter = iter + 1;
    
    disp("-------------------------------");
    disp("Iteration: " + iter-1);
    disp("Initial Deflection Total Error: " + sum(ZError.^2)/ZPointsTotal);
    disp("Initial Pressure Total Error: " + sum(PressError.^2)*ElectAreaPercent);

    ElectInputs = ElectInputs + kappa*PressError;
    %[sortedInputs,idx] = sort(ElectInputs,'descend');    
    sortedInputs = ElectInputs;
    idx = [1:length(ElectInputs)]';    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     ResetFlag = false;
%     for i=1:length(PhiError)
%         if PhiError(i)<0 ResetFlag = true; break; end
%     end
%     
%     if(ResetFlag)
%         disp("Resetting...");
%         title("Resetting");
%         voltageInput = linspace(0,inputMax,samplesPerSegment);
%         [ZVec,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
%         voltageInput = linspace(inputMax,0,samplesPerSegment);
%         [ZVec,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
%         voltageInput = linspace(0,-800,samplesPerSegment);
%         [ZVec,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
%         voltageInput = linspace(-800,0,samplesPerSegment);
%         [ZVec,ZMat,Phi] = applyInputAnim(voltageInput,ElectGrid^2+1,true);
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:length(sortedInputs)
%         title("Iteration: " + (iter+1) + ", Actuating Electrode: " + idx(i));
        voltageInput = [linspace(0,sortedInputs(i),samplesPerSegment)';...
            linspace(sortedInputs(i),0,samplesPerSegment)'];
        [ZVec,ZMat,Press] = applyCoupledInputAnim(voltageInput,idx(i));
        disp("Iteration: " + iter-1 ...
            + ", Step: " + i ...
            + ", Electrode: " + idx(i) ...
            + ", Amplitude: " + sortedInputs(i) ...
            + ", Press: " + Press(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Save figure
%     title("Iteration: " + (iter+1) + ", Actuating Electrode: " + idx(i));
    saveas(gca, "hdm_control_iter_"+iter+".eps","epsc");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

disp("-------------------------------");
disp("Total Iterations: " + iter);
disp("Final Deflection Total Error: " + sum(ZError.^2)/ZPointsTotal);
disp("Final Pressure Total Error: " + sum(PressError.^2)*ElectAreaPercent);

if(isRecording)
    close(vidWriter);
end

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

function [ZVec,ZMat,Press] = applyCoupledInputAnim(input,elect,forceReload)
    persistent PhiArr ECoup HBold MirrorMaskIdxMap MirrorGridSize ...
                surfHand isRecording vidWriter;
%                 ZFit;
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
        % Apply input
        applyInput(PhiArr,input(i),elect,ECoup);
        Press = 1.5*getPhiArrOutputs(PhiArr);

        % Compute deflections
        ZVec = HBold*Press;
        ZMat = MatUtils.vecIdxMapToMatrix(ZVec,MirrorMaskIdxMap,...
                                MirrorGridSize,MirrorGridSize,NaN);
        
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

function applyInput(PhiArr,input,muxCh,ECoup)
    inputs = ECoup(:,muxCh)*input;
    for i=1:length(PhiArr)
        PhiArr(i).updateRelays(inputs(i));
    end
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
