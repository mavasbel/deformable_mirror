% clear all
close all
clc

% Parameters for electrodes creation
lineWidth = 1.5;
MaxRad = 1;
EfectiveSurfFact = 0.70;
ElectGrid = 5;
ElectWidth = 0.12;
ElectCorner00 = EfectiveSurfFact*MaxRad/sqrt(2);
ElectSpace = (EfectiveSurfFact*MaxRad/sqrt(2)-ElectWidth...
                + EfectiveSurfFact*MaxRad/sqrt(2))/(ElectGrid-1);

% Array where each row is an electrode in format (x1 x2 y1 y2)
ElectCorners = zeros(ElectGrid^2,4);
for k = 1:ElectGrid  % These are the rows
    for l = 1:ElectGrid % These are the columns
        ElectCorners(l+(k-1)*ElectGrid,1) = -ElectCorner00            + (l-1)*ElectSpace;
        ElectCorners(l+(k-1)*ElectGrid,2) = -ElectCorner00+ElectWidth + (l-1)*ElectSpace;
        ElectCorners(l+(k-1)*ElectGrid,3) =  ElectCorner00-ElectWidth - (k-1)*ElectSpace;
        ElectCorners(l+(k-1)*ElectGrid,4) =  ElectCorner00            - (k-1)*ElectSpace;
        plot([ElectCorners(l+(k-1)*ElectGrid,1),...
                    ElectCorners(l+(k-1)*ElectGrid,2)],...
                [ElectCorners(l+(k-1)*ElectGrid,3),...
                    ElectCorners(l+(k-1)*ElectGrid,3)],...
                'k','linewidth',lineWidth); hold on;
        plot([ElectCorners(l+(k-1)*ElectGrid,2),...
                    ElectCorners(l+(k-1)*ElectGrid,2)],...
                [ElectCorners(l+(k-1)*ElectGrid,3),...
                    ElectCorners(l+(k-1)*ElectGrid,4)],...
                'k','linewidth',lineWidth); hold on;
        plot([ElectCorners(l+(k-1)*ElectGrid,2),...
                    ElectCorners(l+(k-1)*ElectGrid,1)],...
                [ElectCorners(l+(k-1)*ElectGrid,4),...
                    ElectCorners(l+(k-1)*ElectGrid,4)],...
                'k','linewidth',lineWidth); hold on;
        plot([ElectCorners(l+(k-1)*ElectGrid,1),...
                    ElectCorners(l+(k-1)*ElectGrid,1)],...
                [ElectCorners(l+(k-1)*ElectGrid,4),...
                    ElectCorners(l+(k-1)*ElectGrid,3)],...
                'k','linewidth',lineWidth); hold on;
        if k == 1
            plot([ElectCorners(l+(k-1)*ElectGrid,2),...
                    ElectCorners(l+(k-1)*ElectGrid,2)],...
                [-1,1]*(sqrt( 1-ElectCorners(l+(k-1)*ElectGrid,2)^2 )),...
                '--k','linewidth',0.5*lineWidth); hold on;
            plot([ElectCorners(l+(k-1)*ElectGrid,1),...
                    ElectCorners(l+(k-1)*ElectGrid,1)],...
                [-1,1]*(sqrt( 1-ElectCorners(l+(k-1)*ElectGrid,1)^2 )),...
                '--k','linewidth',0.5*lineWidth); hold on;
        end
    end
    plot([-1,1]*(sqrt( 1-ElectCorners(l+(k-1)*ElectGrid,3)^2 )),...
    [ElectCorners(l+(k-1)*ElectGrid,3),...
        ElectCorners(l+(k-1)*ElectGrid,3)],...
    'k','linewidth',0.5*lineWidth); hold on;
    plot([-1,1]*(sqrt( 1-ElectCorners(l+(k-1)*ElectGrid,4)^2 )),...
    [ElectCorners(l+(k-1)*ElectGrid,4),...
        ElectCorners(l+(k-1)*ElectGrid,4)],...
    'k','linewidth',0.5*lineWidth); hold on;
end
x = linspace(-MaxRad,MaxRad,100);
plot(x,sqrt(MaxRad^2-x.^2),'k','linewidth',lineWidth); hold on;
plot(x,-sqrt(MaxRad^2-x.^2),'k','linewidth',lineWidth); hold on;
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for computation and saving of influence functions
a = 1;
T = 1;
MirrorGridSize = 100;
mirrorXLim = [-MaxRad, MaxRad];
mirrorYLim = [-MaxRad, MaxRad];
mirrorXIsoLine = linspace(mirrorXLim(1),mirrorXLim(2),MirrorGridSize);
mirrorYIsoLine = linspace(mirrorYLim(1),mirrorYLim(2),MirrorGridSize);
[mirrorXGrid,mirrorYGrid] = meshgrid(mirrorXIsoLine,mirrorYIsoLine);

MirrorMask = mirrorXGrid.^2 + mirrorYGrid.^2 <= MaxRad^2*ones(MirrorGridSize,MirrorGridSize);

% for j=1:MirrorGridSize
%     for i=1:MirrorGridSize
%         % Validates that point is inside the mirror
%         if( MirrorMask(i,j) == 1)
%             plot(mirrorXGrid(i,j),mirrorYGrid(i,j),'*b'); hold on;
%             drawnow limitrate
%         end
%     end
% end
