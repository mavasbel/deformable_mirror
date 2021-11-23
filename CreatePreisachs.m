close all
clc

% Consider electrical coupling factor and paramenters of base Preisach
CouplFac = 0.40;
basePhi = preisachRelayModel;
inputMin = preisachRelayModel.inputGrid(1);
inputMax = preisachRelayModel.inputGrid(end);
gridSize = preisachRelayModel.gridSize;

% Creates vector of Preisach operators
PhiArr = [];
for i=1:ElectGrid^2
    Phi = PreisachRelayModel([inputMin, inputMax], gridSize);
    Phi.resetRelaysOff();
    Phi.weightFunc = basePhi.weightFunc;
    Phi.offset = basePhi.offset;
    PhiArr = [PhiArr; Phi];
end

% Creates matrix of coupling factors
ECoup = eye(ElectGrid^2,ElectGrid^2+1);
for i=1:ElectGrid^2
    for j=1:ElectGrid^2
        if(i~=j)
            ECoup(i,j) = CouplFac;
        end
    end
end
ECoup(:,end) = ones(ElectGrid^2,1); % The last corresponds to all selected

% applyInput(PhiArr,inputMax,1,ones(ElectGrid^2,ElectGrid^2));
% applyInput(PhiArr,-800,1,ones(ElectGrid^2,ElectGrid^2));
% applyInput(PhiArr,0,1,ones(ElectGrid^2,ElectGrid^2));
%  [i,getPhiArrOutputs(PhiArr)'/1000]
% 
% steps = 200;
% v = [linspace(0,inputMax,steps)';linspace(inputMax,0,steps)'];
% for i=1:length(v)
%     applyInput(PhiArr,v(i),1,ECoup);
%     [i,getPhiArrOutputs(PhiArr)'/1000]
% end
% 
% function applyInput(PhiArr,input,muxCh,ECoup)
%     inputs = ECoup(muxCh,:)*input;
%     for i=1:length(PhiArr)
%         PhiArr(i).updateRelays(inputs(i));
%     end
% end
% 
% function YPhi = getPhiArrOutputs(PhiArr)
%     YPhi = zeros(length(PhiArr),1);
%     for i=1:length(PhiArr)
%         YPhi(i) = PhiArr(i).getOutput();
%     end
% end