classdef ENVIRONMENT
  properties
    TEMPERATURE;
    HUMIDITY;
    %add new envirronment conditions
  end
  methods
    %----------------------------------------------------------------------
    function this=ENVIRONMENT()
      this.TEMPERATURE =  0.0;
      this.HUMIDITY    =  0.0;
    end
    %----------------------------------------------------------------------
    function this=set.TEMPERATURE(this,value)
      this.TEMPERATURE = value;
    end
    %----------------------------------------------------------------------
    function temp=get.TEMPERATURE(this)
      temp=this.TEMPERATURE;
    end
    %----------------------------------------------------------------------
    function this=set.HUMIDITY(this,value)
      this.HUMIDITY = value;
    end
    %----------------------------------------------------------------------
    function humidity=get.HUMIDITY(this)
      humidity=this.HUMIDITY;
    end
    %----------------------------------------------------------------------
  end
end