close all
clc

% Parameters for computation and saving of influence functions
fileName = "InfMat.mat";
a = 1;
T = 1;
NSeqLim = 80;
ElectGridSize = 50;
MirrorGridSize = 100;
mirrorXLim = [-MaxRad, MaxRad];
mirrorYLim = [-MaxRad, MaxRad];
mirrorXIsoLine = linspace(mirrorXLim(1),mirrorXLim(2),MirrorGridSize);
mirrorYIsoLine = linspace(mirrorYLim(1),mirrorYLim(2),MirrorGridSize);
[mirrorXGrid,mirrorYGrid] = meshgrid(mirrorXIsoLine,mirrorYIsoLine);

tStart = tic;
MirrorMask = mirrorXGrid.^2 + mirrorYGrid.^2 <= MaxRad^2*ones(MirrorGridSize,MirrorGridSize);
InfFuncs = NaN(MirrorGridSize,MirrorGridSize,size(ElectCorners,1));
for k = 1:size(ElectCorners,1)
    tElect = tic;
    for j=1:MirrorGridSize
        parfor i=1:MirrorGridSize
            % Validates that point is inside the mirror
%             if(mirrorXGrid(i,j)^2 +...
%                 mirrorYGrid(i,j)^2 <= MaxRad^2)
            if( MirrorMask(i,j) == 1)
                InfFuncs(i,j,k) = getMTerm(...
                    mirrorXGrid(i,j),...
                    mirrorYGrid(i,j),...
                    ElectCorners(k,1),...
                    ElectCorners(k,2),...
                    ElectCorners(k,3),...
                    ElectCorners(k,4),...
                    ElectGridSize,NSeqLim);
            end
        end
    end
    toc(tElect)
    
    InfFuncs(:,:,k) = (InfFuncs(:,:,k)/ElectGridSize^2)*a^2/(2*pi*T);
    surf(mirrorXGrid,mirrorYGrid,InfFuncs(:,:,k)); hold off;
    title("Electrode:" + k);
    drawnow
end
toc(tStart)

% MCal(:,:,:) = (MCal(:,:,:)/ElectGridSize^2)*a^2/(2*pi*T);

% Save computation and parameters
save(fileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function term = getMTerm(x,y,x1,x2,y1,y2,electGridSize,nlim)
    xIsoLine = linspace(x1,x2,electGridSize);
    yIsoLine = linspace(y1,y2,electGridSize);
    [xGrid,yGrid] = meshgrid(xIsoLine,yIsoLine);
    
    funMat = zeros(electGridSize,electGridSize);
    r = sqrt(x^2+y^2);
    phi = atan2(y,x);
    for j=1:length(xIsoLine)
        for i=1:length(yIsoLine)            
            rp = sqrt(xGrid(i,j)^2 + yGrid(i,j)^2);
            phip = atan2(yGrid(i,j),xGrid(i,j));
            
            %For the change of variables to cartesian coordinates with the
            %transformation map (x,y = T(r,phi) = r*cos(phi),r*sin(phi) if
            %would be neccesary to multiply funMat by 1/det(dT) with
            %1/det([cos(phip) -rp*sin(phip);
            %       sin(phip)  rp*cos(phip)])=1/rp;
            %However the 1/rp cancels with the *rp of the integrand and
            %therefore it is ommited in the code
            if(rp<r)
                funMat(i,j)=getIntegrand1(r,phi,rp,phip,nlim);
            else
                funMat(i,j)=getIntegrand2(r,phi,rp,phip,nlim);
            end
        end
    end
    term = sum(sum(funMat));
end

function int1 = getIntegrand1(r,phi,rp,phip,nlim)
    sum1 = 0;
    for n=1:nlim
        sum1 = sum1 + (1/n)*((r*rp)^n-(rp/r)^n)...
            *cos(n*(phip-phi));
    end
    %Due to the change of variables *rp is not neccesary because it cancels
    %with the 1/det(dT)=1/rp
    int1 = (log(1/r)-sum1); 
end

function int2 = getIntegrand2(r,phi,rp,phip,nlim)
    sum2 = 0;
    for n=1:nlim
        sum2 = sum2 + (1/n)*((r*rp)^n-(r/rp)^n)...
            *cos(n*(phip-phi));
    end
    %Due to the change of variables *rp is not neccesary because it cancels
    %with the 1/det(dT)=1/rp
    int2 = (log(1/rp)-sum2);
end

% function int1 = getIntegrand1(r,phi,rp,phip,nlim)
%     persistent n
%     if isempty(n)
%         syms n;
%     end
%     sum1 = symsum( (1/n)*((r*rp)^n-(rp/r)^n)...
%                     *cos(n*(phip-phi)),...
%                     n,1,nlim );
%     %Due to the change of variables *rp is not neccesary because it cancels
%     %with the 1/det(dT)= 1/rp
%     int1 = (log(1/r)-sum1); 
% end

% function int2 = getIntegrand2(r,phi,rp,phip,nlim)
%     persistent n
%     if isempty(n)
%         syms n;
%     end
%     sum2 = symsum( (1/n)*((r*rp)^n-(r/rp)^n)...
%                     *cos(n*(phip-phi)),...
%                     n,1,nlim );
%     int2 = (log(1/rp)-sum2);
% end