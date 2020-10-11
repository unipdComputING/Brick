classdef GLOBAL
    properties(Constant)
        %costanti
        TOTELNODES = 8;
        INTPOINT   = 8;
        DIMSPACE   = 3;
        TOTSTATEV  = 0;
        NDOF       = 3; %Ux,Uy,Uz
        SIX        = 6;
        %solver parameters
        TOLL       = 0.000001;
        maxIter    = 150 ;
        %material parameters
        TOLLMAT    = 0.00000001;
        MAXITERMAT = 20;
    end
    methods
        %------------------------------------------------------------------
        function obj=GLOBAL()
        end
        %------------------------------------------------------------------
        function [distance]=getDistance(obj,X2,X1)
            distance = 0.0;
            try
                X = X2-X1;
                distance = norm(X);
            catch
                disp('ERROR in cGLOBAL.getDistance: distance not found\n')
            end
        end
        %------------------------------------------------------------------
        
    end
end

