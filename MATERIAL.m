classdef MATERIAL < GLOBAL
  properties
    E;
    nu;
    density;
    alpha;
  end
  methods
    %-----------------------------------------------------------constructor
    function this=MATERIAL(elasticModulus,poissonCoefficient,density,alpha)
      this.E       = 0.0;
      this.nu      = 0.0;
      this.density = 0.0;
      this.alpha   = 0.0;
      if nargin == 1
        this.E       = elasticModulus;
      elseif nargin == 2
        this.E       = elasticModulus;
        this.nu      = poissonCoefficient;
      elseif nargin == 3
        this.E       = elasticModulus;
        this.nu      = poissonCoefficient;
        this.density = density;
      else
        this.E       = elasticModulus;
        this.nu      = poissonCoefficient;
        this.density = density;
        this.alpha   = alpha;
      end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function setElasticModulus(this,value)
      this.E  = value;
    end
    function this=set.E(this,value)
      this.E  = value;
    end
    %----------------------------------------------------------------------
    function [E]=getElasticValue(this)
      E = this.E;
    end
    function E=get.E(this)
      E = this.E;
    end
    %----------------------------------------------------------------------
    function setPoissonCoefficient(this,value)
      this.nu  = value;
    end
    function this=set.nu(this,value)
      this.nu=value;
    end
    %----------------------------------------------------------------------
    function [nu]=getPoissonCoefficient(this)
      nu = this.nu;
    end
    function nu=get.nu(this)
      nu = this.nu;
    end
    %----------------------------------------------------------------------
    function this=set.density(this,value)
      this.density = value;
    end
    %----------------------------------------------------------------------
    function [density]=get.density(this)
      density = this.density;
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function [D]=getMaterialMat(this)
      D  = zeros(this.SIX,this.SIX);
      val01 = this.E/((1.0+this.nu)*(1.0-2.0*this.nu));
      val02 = val01*(1.d0-2.d0*this.nu)/2.0;
      D(1,1) = (1.0-this.nu)*val01; D(1,2) = val01*this.nu;       D(1,3) = val01*this.nu;
      D(2,1) = val01*this.nu;       D(2,2) = (1.0-this.nu)*val01; D(2,3) = val01*this.nu;
      D(3,1) = val01*this.nu;       D(3,2) = val01*this.nu;       D(3,3) = (1.0-this.nu)*val01;
      D(4,4) = val02; 
      D(5,5) = val02; 
      D(6,6) = val02;
    end
    %----------------------------------------------------------------------
    function [D,stress,statev]=getConstitutiveMat(this,dStrain,strain,statev,temp)
      D = this.getMaterialMat();
      eStrain= zeros(this.SIX,1);
      stress = zeros(this.SIX,1);
      T1 = temp(1); T2 = temp(2); dT = T2-T1;
      alpha = this.alpha; %thermal expantion
      %elastic strain
      eStrain(1) = strain(1)+dStrain(1)-alpha*dT;
      eStrain(2) = strain(2)+dStrain(2)-alpha*dT;
      eStrain(3) = strain(3)+dStrain(3)-alpha*dT;
      eStrain(4) = strain(4)+dStrain(4);
      eStrain(5) = strain(5)+dStrain(5);
      eStrain(6) = strain(6)+dStrain(6);
      %Stress
      stress(1)  =D(1,1)*eStrain(1)+D(1,2)*eStrain(2)+D(1,3)*eStrain(3);
      stress(2)  =D(2,1)*eStrain(1)+D(2,2)*eStrain(2)+D(2,3)*eStrain(3);
      stress(3)  =D(3,1)*eStrain(1)+D(3,2)*eStrain(2)+D(3,3)*eStrain(3);
      stress(4)  =D(4,4)*eStrain(4);
      stress(5)  =D(5,5)*eStrain(5);
      stress(6)  =D(6,6)*eStrain(6);
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
  end
end