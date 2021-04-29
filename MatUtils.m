classdef MatUtils
    
    methods(Static)
        function [vec,idxMap] = matrixToVecIdxMap(mat,mask)
            vec = NaN(sum(sum(mask)),1);
            idxMap = NaN(sum(sum(mask)),2);

            idx = 1;
            for j=1:size(mat,2)
                for i=1:size(mat,1)
                    if(mask(i,j)==1)
                        vec( idx ) = mat(i,j);
                        idxMap( idx , : ) = [i,j];
                        idx = idx+1;
                    end
                end
            end
        end
        
        function matrix = vecIdxMapToMatrix(vec,idxMap,rows,cols,fill)
            matrix = ones(rows,cols)*fill;

            for k=1:length(vec)
                matrix(idxMap(k,1),idxMap(k,2)) = vec(k);
            end
        end
    end
    
end