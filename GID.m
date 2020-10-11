classdef GID < GLOBAL
  properties
  end
  methods
    %-----------------------------------------------------------constructor
    function this=GID
    end
    %----------------------------------------------------------------------
    function mshBEAM(this,stFilePath, nodes,elements)
      try
        path = strcat(stFilePath,'.msh');
        unit = fopen(path, 'wt');
        totNodes = size(nodes,1);
        totElem  = size(elements,1);
        fprintf(unit,'#RESULTS Scilab Brick8_01 models \n');
        fprintf(unit,'MESH dimension 3 Elemtype Hexahedra Nnode 8\n');
        fprintf(unit,'#Nodal coordinates \n');
        fprintf(unit,'coordinates \n');
        for n=1:totNodes
            fprintf(unit, '%i %20.12f %20.12f %20.12f \n',n,nodes(n,1), ...
                                                            nodes(n,2), ...
                                                            nodes(n,3));
        end
        %
        fprintf(unit,'end coordinates \n');
        fprintf(unit,'elements \n');
        for e=1:totElem
            fprintf(unit, '%i %i %i %i %i %i %i %i %i \n', e, elements(e).getNodeIndex(1), ...
                                                              elements(e).getNodeIndex(2), ...
                                                              elements(e).getNodeIndex(3), ...
                                                              elements(e).getNodeIndex(4), ...
                                                              elements(e).getNodeIndex(5), ...
                                                              elements(e).getNodeIndex(6), ...
                                                              elements(e).getNodeIndex(7), ...
                                                              elements(e).getNodeIndex(8));
        end
        fprintf(unit,'end elements \n');
        fclose(unit);
      catch ME
        fprintf('WORNING in GID->mshBEAM cannot open file for writing\n');
        st = strcat(ME.identifier,'\n');
        fprintf(st);
      end
    end
    %----------------------------------------------------------------------
    function [unit]=openCloseGIDResFile(this,stFilePath,unit,status)
      try
        if (strcmp(status,'open')==1)
          path = strcat(stFilePath,'.res');
          unit = fopen(path, 'w');
          %
          fprintf(unit,'GiD Post Results File 1.0 \n');
          fprintf(unit,'\n');
          fprintf(unit,'GaussPoints "8PGauss" Elemtype Hexahedra\n');
          fprintf(unit,'Number Of Gauss Points:     8\n');
          fprintf(unit,'Natural Coordinates: given\n');
          fprintf(unit,'-0.577350269189626 -0.577350269189626 -0.577350269189626\n');
          fprintf(unit,' 0.577350269189626 -0.577350269189626 -0.577350269189626\n');
          fprintf(unit,' 0.577350269189626  0.577350269189626 -0.577350269189626\n');
          fprintf(unit,'-0.577350269189626  0.577350269189626 -0.577350269189626\n');
          fprintf(unit,'-0.577350269189626 -0.577350269189626  0.577350269189626\n');
          fprintf(unit,' 0.577350269189626 -0.577350269189626  0.577350269189626\n');
          fprintf(unit,' 0.577350269189626  0.577350269189626  0.577350269189626\n');
          fprintf(unit,'-0.577350269189626  0.577350269189626  0.577350269189626\n');
          fprintf(unit,'end gausspoints\n');
          fprintf(unit,'\n');
        elseif (strcmp(status,'close')==1)
          fclose(unit);
        end
      catch ME
        fprintf('WARNING in GID->openCloseGIDResFile\n');
        st = strcat(ME.identifier,'\n');
        fprintf(st);
      end
    end
    %----------------------------------------------------------------------
    function plotGID_Node_Scalar(this,unit,nameRes,nameAnalisys, ...
                                 nameScalar, ...
                                 time,values)
      try
        gl = GLOBAL();
        st = strcat('Result ',' "',nameRes     ,'" "', ...
                      nameAnalisys,'" ',num2str(time,4),' Scalar OnNodes\n');
        fprintf(unit,st);
        st = strcat('ComponentNames "',nameScalar,'"\n');
        fprintf(unit,st);
        fprintf(unit,'Values\n');
        dim = size(values,1);
        if dim >= 1
          for n=1:dim
            v1 = values(n);
            st = strcat('  ' , num2str(n) , '  %f\n');
            fprintf(unit,st,v1);
          end
        end
        fprintf(unit,'End Values\n');
      catch ME
        fclose(unit);
        st = strcat(ME.identifier,'\n');
        fprintf(st);
      end
    end
    %----------------------------------------------------------------------
    function plotGID_Node_Vector(this,unit,nameRes,nameAnalisys, ...
                             nameDir1,nameDir2,nameDir3, ...
                             time,nodes,vector)
      try
        gl = GLOBAL();
        st = strcat('Result ',' "',nameRes     ,'" "', ...
                      nameAnalisys,'" ',num2str(time,4),' Vector OnNodes\n');
        fprintf(unit,st);
        st = strcat('ComponentNames',' "',nameDir1,'" "',nameDir2,'" "',nameDir3,'"\n');
        fprintf(unit,st);
        fprintf(unit,'Values\n');
        dim = size(nodes,1);
        if dim >= 1
          for n=1:dim
            v1 = vector((n-1)*gl.NDOF+1);
            v2 = vector((n-1)*gl.NDOF+2);
            v3 = vector((n-1)*gl.NDOF+3);
            st = strcat('  ' , num2str(n) , '  %f %f %f\n');
            fprintf(unit,st,v1,v2,v3);
          end
        end
        fprintf(unit,'End Values\n');
      catch ME
        fclose(unit);
        st = strcat(ME.identifier,'\n');
        fprintf(st);
      end
    end
    %----------------------------------------------------------------------
    function plotGID_IP_Matrix(this,unit,nameRes,nameAnalisys, ...
                           nameDir1,nameDir2,nameDir3, ...
                           nameDir4,nameDir5,nameDir6, ...
                           time,elements,matrixRes)
      try
        gl=GLOBAL();
        st = strcat('Result ',' "',nameRes     ,'" "' ...
                      ,nameAnalisys,'" ',' ',num2str(time,4) ...
                      ,' Matrix  OnGaussPoints "8PGauss"\n');
        fprintf(unit,st);
        st = strcat('ComponentNames',' "',nameDir1,'" "',nameDir2,'" "',nameDir3 ...
                             ,'" "',nameDir4,'" "',nameDir5,'" "',nameDir6,'"\n');
        fprintf(unit,st);
        fprintf(unit,'Values\n');
        dim = size(elements,1);
        if dim >= 1 
          for e=1:dim
            v1 = matrixRes((e-1)*gl.INTPOINT+1,1);
            v2 = matrixRes((e-1)*gl.INTPOINT+1,2);
            v3 = matrixRes((e-1)*gl.INTPOINT+1,3);
            v4 = matrixRes((e-1)*gl.INTPOINT+1,4);
            v5 = matrixRes((e-1)*gl.INTPOINT+1,5);
            v6 = matrixRes((e-1)*gl.INTPOINT+1,6);
            %
            st = '  %i  %f %f %f %f %f %f\n';
            fprintf(unit,st,e,v1,v2,v3,v4,v5,v6);
            for ip=2:gl.INTPOINT
              v1 = matrixRes((e-1)*gl.INTPOINT+ip,1);
              v2 = matrixRes((e-1)*gl.INTPOINT+ip,2);
              v3 = matrixRes((e-1)*gl.INTPOINT+ip,3);
              v4 = matrixRes((e-1)*gl.INTPOINT+ip,4);
              v5 = matrixRes((e-1)*gl.INTPOINT+ip,5);
              v6 = matrixRes((e-1)*gl.INTPOINT+ip,6);
              %
              st = '      %f %f %f %f %f %f\n';
              fprintf(unit,st,v1,v2,v3,v4,v5,v6);
            end
          end
        end
        fprintf(unit,'End Values\n');
      catch ME
        fprintf('WARNING: plotGID_IP_Matrix');
        st = strcat(ME.identifier,'\n');
        fprintf(st);
        fclose(unit);
      end
    end
    %----------------------------------------------------------------------
    function plotGiD_STATEV(this,unit,nameAnalisys, ...
                        time,elements,STATEV)
      try
        gl = GLOBAL();
        for nSv=1:gl.TOTSTATEV
          nameRes =  strcat('sv',num2str(nSv));
          st = strcat('Result ',' "',nameRes     ,'" "' ...
                        ,nameAnalisys,'" ',' ',num2str(time,4) ...
                        ,' Scalar  OnGaussPoints "8PGauss"\n');
          fprintf(unit,st);
          nameComp = strcat('sv',num2str(nSv));
          st = strcat('ComponentNames',' "',nameComp,'"\n');
          fprintf(unit,st);
          fprintf(unit,'Values\n');
          dim = size(elements,1);
          if dim >= 1 
            for e=1:dim
              v1 = STATEV((e-1)*gl.INTPOINT+1,nSv);
              st = '  %i  %f\n';
              fprintf(unit,st,e,v1);
              for ip=2:gl.INTPOINT
                v1 = STATEV((e-1)*gl.INTPOINT+ip,nSv);
                st = '      %f\n';
                fprintf(unit,st,v1);
              end
            end
          end
          fprintf(unit,'End Values\n') ; 
        end
      catch ME
        fprintf('WARNING: plotGiD_STATEV');
        st = strcat(ME.identifier,'\n');
        fprintf(st);
        fclose(unit);
      end 
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
  end
end