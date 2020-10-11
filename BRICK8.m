classdef BRICK8 < GLOBAL
  properties
    INTEGRATIONWEIGHTS= 1.0;
    TOTDOFNODES;
    indexMat;
    nodes;
  end
  methods
    %-----------------------------------------------------------constructor
    function this=BRICK8(nodes,indexMat)
      this.TOTDOFNODES= this.NDOF*this.TOTELNODES;
      this.nodes      = zeros(this.TOTELNODES,1);
      this.indexMat   = 0;
      if nargin >0
        this.nodes = nodes;
        this.indexMat = indexMat;
      end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %function setIndexMaterial(this,index)
    %  this.indexMat = index;
    %end
    function this=set.indexMat(this,value)
      this.indexMat=value;
    end
    %----------------------------------------------------------------------
    function [index]=getIndexMaterial(this)
     index = this.indexMat;
    end
    %----------------------------------------------------------------------
    function setNodesIndex(this,nodesIndexVector)
      try
        this.nodes=nodesIndexVector;
      catch
        fprintf('ERROR in BRICK->setNodesIndex:\n');
        fprimtf('nodesIndexVector not compatibile.\n');
      end
    end
    %----------------------------------------------------------------------
    function [nodeIndex]=getNodeIndex(this,index)
      nodeIndex = 0;
      try
        nodeIndex = this.nodes(index);
      catch
      end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function [CSI]=getLocalNode(this,nodeNumber)
      CSI=zeros(this.NDOF,1);
      if ((nodeNumber>this.TOTELNODES)||(nodeNumber<1))
        fprintf('WARNING in BRICK8 -> getLocalNode\n');
        fprintf('node number is not correct\n');
        return;
      end
      switch nodeNumber
      case 1 
        CSI(1) = -1.0;    CSI(2) = -1.0;    CSI(3) = -1.0;
      case 2 
        CSI(1) =  1.0;    CSI(2) = -1.0;    CSI(3) = -1.0;
      case 3 
        CSI(1) =  1.0;    CSI(2) =  1.0;    CSI(3) = -1.0;
      case 4
        CSI(1) = -1.0;    CSI(2) =  1.0;    CSI(3) = -1.0;
      case 5 
        CSI(1) = -1.0;    CSI(2) = -1.0;    CSI(3) =  1.0;
      case 6 
        CSI(1) =  1.0;    CSI(2) = -1.0;    CSI(3) =  1.0;
      case 7 
        CSI(1) =  1.0;    CSI(2) =  1.0;    CSI(3) =  1.0;
      case 8 
        CSI(1) = -1.0;    CSI(2) =  1.0;    CSI(3) =  1.0;
      end
    end
    %----------------------------------------------------------------------
    function [sf]=getShapeFunction(this, sfNumber, point)
      sf = 0.0;
      if ((sfNumber>this.TOTELNODES)||(sfNumber<1))
        fprintf('WARNING in BRICK8 -> getShapeFunction\n');
        return;
      end
      switch(sfNumber)
      case (1) 
        sf =(1.0-point(1))*(1.0-point(2))*(1.0-point(3))/8.0;
      case (2) 
        sf =(1.0+point(1))*(1.0-point(2))*(1.0-point(3))/8.0;
      case (3) 
        sf =(1.0+point(1))*(1.0+point(2))*(1.0-point(3))/8.0;
      case (4) 
        sf =(1.0-point(1))*(1.0+point(2))*(1.0-point(3))/8.0;
      case (5) 
        sf =(1.0-point(1))*(1.0-point(2))*(1.0+point(3))/8.0;
      case (6) 
        sf =(1.0+point(1))*(1.0-point(2))*(1.0+point(3))/8.0;
      case (7) 
        sf =(1.0+point(1))*(1.0+point(2))*(1.0+point(3))/8.0;
      case (8) 
        sf =(1.0-point(1))*(1.0+point(2))*(1.0+point(3))/8.0;
      end
    end
    %----------------------------------------------------------------------
    function [ipoint]=getIntegrationPoint(this,intPoinNumber)
      ipoint=zeros(3,1);
      if this.INTPOINT==1
          ipoint(1)=  0.000000000000000;
          ipoint(2)=  0.000000000000000;
          ipoint(3)=  0.000000000000000;
      else
        switch  intPoinNumber
        case(1) 
          ipoint(1)= -0.577350269189626;
          ipoint(2)= -0.577350269189626;
          ipoint(3)= -0.577350269189626;
        case(2) 
          ipoint(1)=  0.577350269189626;
          ipoint(2)= -0.577350269189626;
          ipoint(3)= -0.577350269189626;
        case(3) 
          ipoint(1)=  0.577350269189626;
          ipoint(2)=  0.577350269189626;
          ipoint(3)= -0.577350269189626;
        case(4) 
          ipoint(1)= -0.577350269189626;
          ipoint(2)=  0.577350269189626;
          ipoint(3)= -0.577350269189626;
        case(5) 
          ipoint(1)= -0.577350269189626;
          ipoint(2)= -0.577350269189626;
          ipoint(3)=  0.577350269189626;
        case(6) 
          ipoint(1)=  0.577350269189626;
          ipoint(2)= -0.577350269189626;
          ipoint(3)=  0.577350269189626;
        case(7) 
          ipoint(1)=  0.577350269189626;
          ipoint(2)=  0.577350269189626;
          ipoint(3)=  0.577350269189626;
        case(8) 
          ipoint(1)= -0.577350269189626;
          ipoint(2)=  0.577350269189626;
          ipoint(3)=  0.577350269189626;
        end
      end
    end
    %----------------------------------------------------------------------
    function [dSf]=getDerShapeFunction(this,sfNumber, direction, point)
      dSf = 0.0;
      if ((sfNumber>this.TOTELNODES)||(sfNumber<1)) 
        fprintf('WARNING in getDerShapeFunction\n');
        return;
      end
      if ((direction>3)||(direction<1)) 
        fprintf('WARNING in getDerShapeFunction\n');
        return;
      end
      %
      %localNode = this.getLocalNode(sfNumber);
      %
      switch (sfNumber)
      case (1)  %sf 1
        switch (direction)
        case (1) 
          dSf = -1.0*(1.0-point(2))*(1.0-point(3))/8.0;
        case (2) 
          dSf = -1.0*(1.0-point(1))*(1.0-point(3))/8.0;
        case (3) 
          dSf = -1.0*(1.0-point(1))*(1.0-point(2))/8.0;                    
        end
      case (2) %sf 2
        switch (direction)
        case (1) 
          dSf =   1.0*(1.0-point(2))*(1.0-point(3))/8.0;
        case (2) 
          dSf =  -1.0*(1.0+point(1))*(1.0-point(3))/8.0;
        case (3) 
          dSf =  -1.0*(1.0+point(1))*(1.0-point(2))/8.0;
        end
      case (3) %sf 3
        switch (direction)
        case (1) 
          dSf =   1.0*(1.0+point(2))*(1.0-point(3))/8.0;
        case (2) 
          dSf =   1.0*(1.0+point(1))*(1.0-point(3))/8.0;
        case (3) 
          dSf =  -1.0*(1.0+point(1))*(1.0+point(2))/8.0;
        end
      case (4)  %sf 4
        switch (direction)
        case (1) 
          dSf =  -1.0*(1.0+point(2))*(1.0-point(3))/8.0;
        case (2) 
          dSf =   1.0*(1.0-point(1))*(1.0-point(3))/8.0;
        case (3) 
          dSf =  -1.0*(1.0-point(1))*(1.0+point(2))/8.0;
        end
      case (5)  %sf 5
        switch (direction)
        case (1) 
          dSf = -1.0*(1.0-point(2))*(1.0+point(3))/8.0;
        case (2) 
          dSf = -1.0*(1.0-point(1))*(1.0+point(3))/8.0;
        case (3) 
          dSf =  1.0*(1.0-point(1))*(1.0-point(2))/8.0;                    
        end
      case (6)  %sf 6
        switch (direction)
        case (1) 
          dSf =   1.0*(1.0-point(2))*(1.0+point(3))/8.0;
        case (2) 
          dSf =  -1.0*(1.0+point(1))*(1.0+point(3))/8.0;
        case (3) 
          dSf =   1.0*(1.0+point(1))*(1.0-point(2))/8.0;
        end
      case (7) %sf 7
        switch (direction)
        case (1) 
          dSf =   1.0*(1.0+point(2))*(1.0+point(3))/8.0;
        case (2) 
          dSf =   1.0*(1.0+point(1))*(1.0+point(3))/8.0;
        case (3) 
          dSf =   1.0*(1.0+point(1))*(1.0+point(2))/8.0;
        end
      case (8) %sf 8
        switch (direction)
        case (1) 
          dSf =  -1.0*(1.0+point(2))*(1.0+point(3))/8.0;
        case (2) 
          dSf =   1.0*(1.0-point(1))*(1.0+point(3))/8.0;
        case (3) 
          dSf =   1.0*(1.0-point(1))*(1.0+point(2))/8.0;
        end
      end
    end
    %----------------------------------------------------------------------
    function [dSf]=getDerShapeFunction01(this,direction, point)
      dSf = zeros(8,1);
      if ((direction>3)||(direction<1))
        fprintf('WARNING in getDerShapeFunction\n');
        return;
      end
      %
      %localNode = getLocalNode(sfNumber);
      %
      switch (direction)
      case (1) 
        dSf(1) =  -1.0*(1.0-point(2))*(1.0-point(3))/8.0;
        dSf(2) =   1.0*(1.0-point(2))*(1.0-point(3))/8.0;
        dSf(3) =   1.0*(1.0+point(2))*(1.0-point(3))/8.0;
        dSf(4) =  -1.0*(1.0+point(2))*(1.0-point(3))/8.0;
        dSf(5) =  -1.0*(1.0-point(2))*(1.0+point(3))/8.0;
        dSf(6) =   1.0*(1.0-point(2))*(1.0+point(3))/8.0;
        dSf(7) =   1.0*(1.0+point(2))*(1.0+point(3))/8.0;
        dSf(8) =  -1.0*(1.0+point(2))*(1.0+point(3))/8.0;
      case (2) 
        dSf(1) =  -1.0*(1.0-point(1))*(1.0-point(3))/8.0;
        dSf(2) =  -1.0*(1.0+point(1))*(1.0-point(3))/8.0;
        dSf(3) =   1.0*(1.0+point(1))*(1.0-point(3))/8.0;
        dSf(4) =   1.0*(1.0-point(1))*(1.0-point(3))/8.0;
        dSf(5) =  -1.0*(1.0-point(1))*(1.0+point(3))/8.0;
        dSf(6) =  -1.0*(1.0+point(1))*(1.0+point(3))/8.0;
        dSf(7) =   1.0*(1.0+point(1))*(1.0+point(3))/8.0;
        dSf(8) =   1.0*(1.0-point(1))*(1.0+point(3))/8.0;
      case (3)  
        dSf(1) =  -1.0*(1.0-point(1))*(1.0-point(2))/8.0;
        dSf(2) =  -1.0*(1.0+point(1))*(1.0-point(2))/8.0;
        dSf(3) =  -1.0*(1.0+point(1))*(1.0+point(2))/8.0;
        dSf(4) =  -1.0*(1.0-point(1))*(1.0+point(2))/8.0;
        dSf(5) =   1.0*(1.0-point(1))*(1.0-point(2))/8.0;
        dSf(6) =   1.0*(1.0+point(1))*(1.0-point(2))/8.0;
        dSf(7) =   1.0*(1.0+point(1))*(1.0+point(2))/8.0;
        dSf(8) =   1.0*(1.0-point(1))*(1.0+point(2))/8.0;
      end
    end
    %----------------------------------------------------------------------
    function [matB]=getMatrixB(this,sfNumber,localPoint,invJacob)
      matB=zeros(this.SIX,this.DIMSPACE);
      dSf1 = this.getDerShapeFunction(sfNumber,1,localPoint);
      dSf2 = this.getDerShapeFunction(sfNumber,2,localPoint);
      dSf3 = this.getDerShapeFunction(sfNumber,3,localPoint);
      dN(1) = invJacob(1,1)*dSf1+...
              invJacob(1,2)*dSf2+...
              invJacob(1,3)*dSf3;
      dN(2) = invJacob(2,1)*dSf1+...
              invJacob(2,2)*dSf2+...
              invJacob(2,3)*dSf3;
      dN(3) = invJacob(3,1)*dSf1+...
              invJacob(3,2)*dSf2+...
              invJacob(3,3)*dSf3;       
      matB(1,1) = dN(1);
      matB(2,2) = dN(2);
      matB(3,3) = dN(3);
      matB(4,1) = dN(2);
      matB(4,2) = dN(1);
      matB(5,2) = dN(3);
      matB(5,3) = dN(2);
      matB(6,1) = dN(3);
      matB(6,3) = dN(1);
    end
    %----------------------------------------------------------------------
    function [invJacob,vol]=jacobian(this,X,localPoint)
      deter    = 0.0;
      vol      = 0.0;
      dN      =zeros(this.DIMSPACE  ,this.TOTELNODES);
      jacob   =zeros(this.DIMSPACE  ,this.DIMSPACE);
      invJacob=zeros(this.DIMSPACE  ,this.DIMSPACE);
      %
      dN(1,:) = this.getDerShapeFunction01(1, localPoint);
      %
      dN(2,:) = this.getDerShapeFunction01(2, localPoint);
      %
      dN(3,:) = this.getDerShapeFunction01(3, localPoint);
      %
      jacob(1,1)  = dN(1,1:this.TOTELNODES)*X(1:this.TOTELNODES,1);
      jacob(2,1)  = dN(2,1:this.TOTELNODES)*X(1:this.TOTELNODES,1);
      jacob(3,1)  = dN(3,1:this.TOTELNODES)*X(1:this.TOTELNODES,1);
      jacob(1,2)  = dN(1,1:this.TOTELNODES)*X(1:this.TOTELNODES,2);
      jacob(2,2)  = dN(2,1:this.TOTELNODES)*X(1:this.TOTELNODES,2);
      jacob(3,2)  = dN(3,1:this.TOTELNODES)*X(1:this.TOTELNODES,2);
      jacob(1,3)  = dN(1,1:this.TOTELNODES)*X(1:this.TOTELNODES,3);
      jacob(2,3)  = dN(2,1:this.TOTELNODES)*X(1:this.TOTELNODES,3);
      jacob(3,3)  = dN(3,1:this.TOTELNODES)*X(1:this.TOTELNODES,3);

      deter = +jacob(1,1)*(jacob(2,2)*jacob(3,3)-jacob(2,3)*jacob(3,2))...
              -jacob(1,2)*(jacob(2,1)*jacob(3,3)-jacob(2,3)*jacob(3,1))...
              +jacob(1,3)*(jacob(2,1)*jacob(3,2)-jacob(2,2)*jacob(3,1));

      if (deter <= 0 )
          %fprintf('ERROR in Jacobian function: determinant is <0\n');
          error('ERROR in Jacobian function: determinant is <0\n');
      end
      %
      %Jacobian^-1 matrix $\textbf{J}^{-1}$
      invJacob(1,1)= 1/deter*(jacob(2,2)*jacob(3,3)-jacob(2,3)*jacob(3,2));
      invJacob(2,1)=-1/deter*(jacob(2,1)*jacob(3,3)-jacob(2,3)*jacob(3,1));
      invJacob(3,1)= 1/deter*(jacob(2,1)*jacob(3,2)-jacob(2,2)*jacob(3,1));
      invJacob(1,2)=-1/deter*(jacob(1,2)*jacob(3,3)-jacob(1,3)*jacob(3,2));
      invJacob(2,2)= 1/deter*(jacob(1,1)*jacob(3,3)-jacob(1,3)*jacob(3,1));
      invJacob(3,2)=-1/deter*(jacob(1,1)*jacob(3,2)-jacob(1,2)*jacob(3,1));
      invJacob(1,3)= 1/deter*(jacob(1,2)*jacob(2,3)-jacob(1,3)*jacob(2,2));
      invJacob(2,3)=-1/deter*(jacob(1,1)*jacob(2,3)-jacob(1,3)*jacob(2,1));
      invJacob(3,3)= 1/deter*(jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1));
      %
      vol = deter*this.INTEGRATIONWEIGHTS;
    end
    %----------------------------------------------------------------------
    function [K,elFi,elStress,STATEV]=localStiffness(this,mat, ...
                                           X,u,DSTRAIN,STRAIN, ...
                                            STATEV,T,time)
      D       =zeros(this.SIX                 ,        this.SIX);
      stress  =zeros(this.SIX                 ,               1);
      K       =zeros(this.TOTDOFNODES         ,this.TOTDOFNODES);
      elFi    =zeros(this.TOTELNODES*this.NDOF,               1);
      elStress=zeros(this.INTPOINT            ,        this.SIX);
      strain  =zeros(this.SIX                 ,               1);
      dStrain =zeros(this.SIX                 ,               1);
      statev  =zeros(this.TOTSTATEV           ,               1);
      BD      =zeros(this.DIMSPACE            , 2*this.DIMSPACE);
      kij     =zeros(this.DIMSPACE            ,   this.DIMSPACE);
      %
      elTemp = this.getElemTemperature(T);
      %
      %integration point loop
      for g=1:this.INTPOINT
        strain     =  STRAIN(g,:);
        dStrain    = DSTRAIN(g,:);
        statev     = STATEV (g,:);
        [D,stress,statev] = mat.getConstitutiveMat(dStrain,strain,statev,elTemp(g,:));
        elStress(g,:)=stress;
        STATEV  (g,:)=statev;
        %
        gaussPoint = this.getIntegrationPoint(g);
        %
        [invJacob,volElement]=this.jacobian(X,gaussPoint);
        %
        for i=1:this.TOTELNODES
          matBi = this.getMatrixB(i,gaussPoint,invJacob);
          for j=1:this.TOTELNODES
            matBj = this.getMatrixB(j,gaussPoint,invJacob);
            %
            BD  = matBi'*D;
            kij = BD*matBj;
            %
            K(3*(i-1)+1,3*(j-1)+1)=K(3*(i-1)+1,3*(j-1)+1)+kij(1,1)*volElement;
            K(3*(i-1)+1,3*(j-1)+2)=K(3*(i-1)+1,3*(j-1)+2)+kij(1,2)*volElement;
            K(3*(i-1)+1,3*(j-1)+3)=K(3*(i-1)+1,3*(j-1)+3)+kij(1,3)*volElement;

            K(3*(i-1)+2,3*(j-1)+1)=K(3*(i-1)+2,3*(j-1)+1)+kij(2,1)*volElement;
            K(3*(i-1)+2,3*(j-1)+2)=K(3*(i-1)+2,3*(j-1)+2)+kij(2,2)*volElement;
            K(3*(i-1)+2,3*(j-1)+3)=K(3*(i-1)+2,3*(j-1)+3)+kij(2,3)*volElement;

            K(3*(i-1)+3,3*(j-1)+1)=K(3*(i-1)+3,3*(j-1)+1)+kij(3,1)*volElement;
            K(3*(i-1)+3,3*(j-1)+2)=K(3*(i-1)+3,3*(j-1)+2)+kij(3,2)*volElement;
            K(3*(i-1)+3,3*(j-1)+3)=K(3*(i-1)+3,3*(j-1)+3)+kij(3,3)*volElement;
            %
          end
          for dof=1:this.NDOF
            matBi = this.getMatrixB(i,gaussPoint,invJacob);
            elFi((i-1)*this.NDOF+dof)=elFi((i-1)*this.NDOF+dof)     + ...
                                   matBi(1,dof)*stress(1)*volElement+ ...
                                   matBi(2,dof)*stress(2)*volElement+ ...
                                   matBi(3,dof)*stress(3)*volElement+ ...
                                   matBi(4,dof)*stress(4)*volElement+ ...
                                   matBi(5,dof)*stress(5)*volElement+ ...
                                   matBi(6,dof)*stress(6)*volElement;
          end
        end
      end
    end
    %----------------------------------------------------------------------
    function [strain]=getElStrain(this,X,nodeDisp)
      %
      strain=zeros(this.INTPOINT,this.SIX);
      for g=1:this.INTPOINT
        gaussPoint     = this.getIntegrationPoint(g);
        [invJacob,vol] = this.jacobian(X,gaussPoint);

        for i=1:this.TOTELNODES
          matB = this.getMatrixB(i,gaussPoint,invJacob);
          strain(g,1) = strain(g,1)+matB(1,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(1,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(1,3)*nodeDisp((i-1)*3+3);
          strain(g,2) = strain(g,2)+matB(2,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(2,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(2,3)*nodeDisp((i-1)*3+3);
          strain(g,3) = strain(g,3)+matB(3,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(3,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(3,3)*nodeDisp((i-1)*3+3);
          strain(g,4) = strain(g,4)+matB(4,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(4,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(4,3)*nodeDisp((i-1)*3+3);
          strain(g,5) = strain(g,5)+matB(5,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(5,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(5,3)*nodeDisp((i-1)*3+3);
          strain(g,6) = strain(g,6)+matB(6,1)*nodeDisp((i-1)*3+1) ...
                                   +matB(6,2)*nodeDisp((i-1)*3+2) ...
                                   +matB(6,3)*nodeDisp((i-1)*3+3);
        end
      end
      %
    end
    %----------------------------------------------------------------------
    function [intElLoad]=getInternalLoad(this,X,stress)
      intElLoad=zeros(this.TOTELNODES*this.NDOF,1);
      for iG=1:this.INTPOINT
        gPoint = this.getIntegrationPoint(iG);
        for n=1:this.TOTELNODES
          [invJacob,vol]=this.jacobian(X,gPoint);
          B =  this.getMatrixB(n,gPoint,invJacob);
          for dof=1:this.NDOF
            intElLoad((n-1)*this.NDOF+dof)=...
                                        intElLoad((n-1)*this.NDOF+dof)+ ...
                                             B(1,dof)*stress(iG,1)*vol+ ...
                                             B(2,dof)*stress(iG,2)*vol+ ...
                                             B(3,dof)*stress(iG,3)*vol+ ...
                                             B(4,dof)*stress(iG,4)*vol+ ...
                                             B(5,dof)*stress(iG,5)*vol+ ...
                                             B(6,dof)*stress(iG,6)*vol;
          end
        end
      end
    end
    %----------------------------------------------------------------------
    function [elTemp]=getElemTemperature(this,nodeTemp)
      INTPOINT  =this.INTPOINT;
      TOTELNODES=this.TOTELNODES;
      elTemp    =zeros(INTPOINT,2);
      %***************************************SELECTIVE REDUCED INTEGRATION
      %x prevenire fenomeni di looking la termica viene calcolata su un
      %punto gauss (si dovrebbe fare anche per la parte idrostatica del
      %tensore delle tensioni)
      point = [0.0;0.0;0.0;];
      for n=1:TOTELNODES
        for ip=1:INTPOINT 
          sf = this.getShapeFunction(n, point);
          elTemp(ip,1) = elTemp(ip,1)+sf*nodeTemp(n,1);
          elTemp(ip,2) = elTemp(ip,2)+sf*nodeTemp(n,2);
        end
      end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
  end
end