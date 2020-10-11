function [nodes,...
          elements,...
          materials,...
          boundaries,...
          loads,...
          temperatures,...
          tables,...
          env,...
          solverPar,...
          filePath]=READ_INPUT_FILE()
  [stPath,pathName] = uigetfile('*.inp');
  stPath = strcat(pathName,stPath);
  try
    unitIN = fopen(stPath,'r');
  catch
    fprintf('ERROR in READ_INPUT_FILE: cannot open input file\n');
    return
  end
  bc = 0;
  m  = 0;
  l  = 0;
  n  = 0;
  e  = 0;
  t  = 0;
  te = 0;
  filePath       =stPath;
  nodes          =[];
  elements(1,1)  =BRICK8([],0);
  materials(1,1) =MATERIAL(0.0,0.0);
  loads          =[];  
  temperatures   =[];
  tables(1,1)    =TABLE(1,[]);
  env            =ENVIRONMENT();
  while ~feof(unitIN)
    stCommand = fgetl(unitIN);
    %
    switch stCommand
    case 'materials'
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        m  = m+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        %
        dim= size(st,2);
        if     dim==2
          E =str2double(st(2)); %E
          materials(m,1)=MATERIAL(E);
        elseif dim==3
          E =str2double(st(2)); %E
          nu=str2double(st(3)); %poisson
          materials(m,1)=MATERIAL(E,nu);
        elseif dim==4
          E    =str2double(st(2)); %E
          nu   =str2double(st(3)); %poisson
          rho  =str2double(st(4)); %density
          materials(m,1)=MATERIAL(E,nu,rho);
        else
          E    =str2double(st(2)); %E
          nu   =str2double(st(3)); %poisson
          rho  =str2double(st(4)); %density
          alpha=str2double(st(5)); %alpha
          materials(m,1)=MATERIAL(E,nu,rho,alpha);
        end
        stLine = fgetl(unitIN);
      end
    case 'boundaries'
      %
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        bc  = bc+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        boundaries(bc,1)=str2double(st(1)); 
        boundaries(bc,2)=str2double(st(2)); 
        boundaries(bc,3)=str2double(st(3));
        boundaries(bc,4)=str2double(st(4)); 
        boundaries(bc,5)=str2double(st(5)); 
        boundaries(bc,6)=str2double(st(6));
        boundaries(bc,7)=str2double(st(7));
        %table
        if size(st,2) >= 8
          boundaries(bc,8)=str2double(st(8));
        end
        stLine = fgetl(unitIN);
      end
    case 'loads'
      %fprintf('loads\n')
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        l  = l+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        loads(l,1)=str2double(st(1)); 
        loads(l,2)=str2double(st(2)); 
        loads(l,3)=str2double(st(3));
        loads(l,4)=str2double(st(4));
        %table
        if (size(st,2)>=5) 
          loads(l,5)=str2double(st(5));
        end
        stLine = fgetl(unitIN);
      end
    case 'temperatures'
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        te  = te+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        temperatures(te,1)=str2double(st(1)); 
        temperatures(te,2)=str2double(st(2)); 
        %table
        if (size(st,2)>=3) 
          temperatures(te,3)=str2double(st(3));
        end
        stLine = fgetl(unitIN);
      end
    case 'nodes'
      %fprintf('nodes\n');
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        n  = n+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = []; 
        offset = 0;
        if (size(st,2) > 3)
          offset = 1;
        end
        nodes(n,1)=str2double(st(1+offset));
        nodes(n,2)=str2double(st(2+offset));
        nodes(n,3)=str2double(st(3+offset));
        %
        stLine = fgetl(unitIN);
      end
    case 'elements'
      %fprintf('elements\n');
      stLine = fgetl(unitIN);
      %
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        e  = e +1;
        indexVal = zeros(9,1);
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        offset = 0;
        if (size(st,2)==10)
          offset = 1;
        end
        indexVal(1)=str2double(st(1+offset));
        indexVal(2)=str2double(st(2+offset));
        indexVal(3)=str2double(st(3+offset));
        indexVal(4)=str2double(st(4+offset));
        indexVal(5)=str2double(st(5+offset));
        indexVal(6)=str2double(st(6+offset));
        indexVal(7)=str2double(st(7+offset));
        indexVal(8)=str2double(st(8+offset));
        indexVal(9)=str2double(st(9+offset));
        elements(e,1)=BRICK8(indexVal(1:8),indexVal(9));
        stLine = fgetl(unitIN);
      end
    case 'table'
      %fprintf('table\n');
      t = t+1;
      tables(t,1) = TABLE(0,[]);
      stLine = fgetl(unitIN);
      if (feof(unitIN)==0)&&(strcmp(stLine,'end')==0)
        %tables(t,1).index(stLine);
        tName = stLine;
        i      = 0;
        stLine = fgetl(unitIN);
        valXY=[];
        while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
          i = i+1;
          st = strsplit(stLine);
          st(strcmp('',st)) = [];
          valXY(i,1) = str2double(st(1));
          valXY(i,2) = str2double(st(2));
          stLine = fgetl(unitIN);
        end
        tables(t,1)= TABLE(tName,valXY);
      end
    case 'solver'
      %fprintf('solverparameters\n');
      i=0;
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        i  = i+1;
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        solverPar(1)=str2double(st(1));
        solverPar(2)=str2double(st(2));
        solverPar(3)=str2double(st(3));
        solverPar(4)=str2double(st(4));
        stLine = fgetl(unitIN);
      end
      %mseek(-10,unitIN,'cur')
    case 'environment'
      stLine = fgetl(unitIN);
      while (~feof(unitIN))&&(strcmp(stLine,'end')==0)
        st = strsplit(stLine);
        st(strcmp('',st)) = [];
        try
          if strcmp(st(1,1),'temp')==1
            env.TEMPERATURE=str2double(st(1,2));
          elseif strcmp(st(1,1),'hum')==1
            env.HUMIDITY=str2double(st(1,2));
          end
        catch
        end
        stLine = fgetl(unitIN);
      end
    end
  end
  fclose(unitIN);
end