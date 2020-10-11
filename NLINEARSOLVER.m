%--------------------------------------------------------------------------
function NLINEARSOLVER(unit,timeStart,timeEnd,totStep,saveEvery,...
                         elements, nodes, materials,            ...
                         loads, temperature, boundary, tables, env )
  %import GLOBAL.*;
  gl      = GLOBAL();
  Fe      = [];
  res     = [];
  strain  = [];
  u       = [];
  Fi      = [];
  du      = [];
  K       = [];
  %
  if nargin>0
    if  totStep < 1
      totStep = 1;
    end
    dimNode     = size(nodes,1);
    dimElem     = size(elements,1);
    dimProblem  = dimNode*gl.NDOF;
    dimIP       = dimElem*gl.INTPOINT;
    dTime       = (timeEnd-timeStart)/totStep;
    time        =  timeStart;
    totStepSaved= fix(totStep/saveEvery);
    %//
    %//
    %//inizializzazione variabili di stato (modelli costitutivi non lin)
    STATEV= zeros(dimElem*dimIP,gl.TOTSTATEV);
    %//
    time= 0.0;
    TIME= zeros(2         ,1);
    u   = zeros(dimProblem,1);
    du  = zeros(dimProblem,1);
    vN  = zeros(dimProblem,1);
    Fe  = zeros(dimProblem,1);
    res = zeros(dimProblem,1);
    %//
     STRESS=zeros(dimElem*gl.INTPOINT,gl.SIX);
     STRAIN=zeros(dimElem*gl.INTPOINT,gl.SIX);
    DSTRAIN=zeros(dimElem*gl.INTPOINT,gl.SIX);
    %
    fprintf('NON LINEAR STATIC SOLVER\n');
    iStepSaved = 0;
    stepSaved  = 1;
    flag = 0;
    %
    %
    %
    vBC = getBCvector(boundary,dimProblem);
    vT  = getTemperatureVector(nodes,temperature,env);
    %
    %
    %
    TIME(2)=timeStart;
    for iStep=1:totStep+1
      fprintf('--Time %10.4f--\n',time);
      solverStatus = 'new';
      Norm = 1.0+gl.TOLL;
      iter = 0;
      %
      du      = add_du_from_BC(du,boundary,tables,TIME);
      %
      DSTRAIN = updateSTRAIN(nodes,elements,du);
      %
      %---
      %external load
      Fe = loadAssembly(nodes,loads,tables,TIME);
      %
      %
      %temperature
      Temp = tempAssembly(vT,env,tables,TIME);
      %
      while Norm>gl.TOLL
        iter=iter+1;
        %assembly (K and Fi)
        [K,Fi,STRESS,STATEV] = assembly(nodes,elements,materials,  ...
                                                  u,STRAIN,DSTRAIN,...
                                                  STATEV,Temp,time,...
                                                  solverStatus);
        %
        res  = (Fe-Fi);
        %
        Norm = getRedidueNorm(res,vBC,iter);
        %
        %solver
        ddu  = solver(K,res,vBC,boundary);
        %aggiornamento
        du  = du+ddu;
        %
        DSTRAIN = updateSTRAIN(nodes,elements,du);
        %
        solverStatus = 'old';
        %
        if iter>gl.maxIter
          fprintf('WARNING: in nonLinearStatic: raggiunto il numero massimo di iterazioni\n');
          fprintf('( |Res|: %10.4e,)\n',Norm);
          fprintf('scrivi resume per terminare il programma\n');
          pause
          fprintf('convergenza non trovata END PROGRAM\n');
          flag=1;
          break
        end
      end
      %
      u      = u+du;
      STRAIN = STRAIN+DSTRAIN;
      %
       if iStep==stepSaved
         iStepSaved=iStepSaved+1;
         gid = GID();
         gid.plotGID_Node_Vector(unit,'Displacements','Static', ...
                                        'X','Y','Z', time,        ...
                                                      nodes,u    );
         gid.plotGID_Node_Scalar(unit,'Temperature','Static', ...
                                               'T',time,Temp(:,2));
         gid.plotGID_IP_Matrix(unit,'Strain','Static', ...
                                 'e11','e22','e33', ...
                                 'e12','e23','e13', ...
                                 time,elements,STRAIN);
         gid.plotGID_IP_Matrix(unit,'Stress','Static', ...
                                 's11','s22','s33', ...
                                 's12','s23','s13', ...
                                 time,elements,STRESS);
         gid.plotGiD_STATEV    (unit,'Static',         ...
                                 time,elements,STATEV);
         stepSaved = stepSaved+saveEvery;
       end
      fprintf(' tot Iter: %5i \n',iter);
      fprintf(' |Res|: %10.4e\n',Norm);
      fprintf('------------------\n');
      if flag==1
        return
      end
      du = zeros(dimProblem,1);
      res= zeros(dimProblem,1);
      time=time+dTime;
      TIME(1)=TIME(2);
      TIME(2)=TIME(2)+dTime;
    end
  end
end
%--------------------------------------------------------------------------
%assemblaggio matrice di rigidezza
%K matrice assemblata
%elK matrice di rigidezza dell'elemento rispetto il s.d.r. globale
function [K,Fi,STRESS,STATEV]=assembly(nodes,elements,                ...
                                       materials,                     ...
                                       displacements,STRAIN, DSTRAIN, ...
                                       STATEV,Temp,time,solverStatus)
  gl= GLOBAL();
  nElnode = gl.TOTELNODES; %numero di nodi per elemento
  dimEl = size(elements,1);
  if dimEl==0
      return;
  end
  dimNode = size(nodes,1);
  if dimNode==0
      return;
  end
  %
  K      =zeros(dimNode*gl.NDOF    ,dimNode*gl.NDOF);
  Fi     =zeros(dimNode*gl.NDOF    ,              1);
  STRESS =zeros(dimEl  *gl.INTPOINT,         gl.SIX);
  dStrain=zeros(gl.INTPOINT        ,         gl.SIX);
   strain=zeros(gl.INTPOINT        ,         gl.SIX);
   statev=zeros(gl.INTPOINT        ,   gl.TOTSTATEV);
  %
  X=zeros(gl.TOTELNODES,gl.NDOF);
  u=zeros(gl.TOTELNODES,gl.NDOF);
  T=zeros(gl.TOTELNODES,      2);
  %
  for e=1:dimEl
    %
    for n=1:gl.TOTELNODES
      indexNode  = elements(e).getNodeIndex(n);
      X(n,:)     = nodes(indexNode,:);
      T(n,:)     = Temp (indexNode,:);
      for dof=1:gl.NDOF
        u(n,dof)= displacements((indexNode-1)*gl.NDOF+dof)';
      end
    end
    m    = elements (e).getIndexMaterial();%indice del materiale
    mat  = materials(m);
    %matrice di rigidezza dell'elemento
      strain((1:gl.INTPOINT),:)= STRAIN((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:);
     dStrain((1:gl.INTPOINT),:)=DSTRAIN((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:);
     statev ((1:gl.INTPOINT),:)= STATEV((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:);
    %
    %
    [elK,elFi,elStress,statev] = elements(e).localStiffness(mat, ...
                                             X,u,dStrain,strain, ...
                                             statev,T,time);
                                         
    %stress assembly
    STRESS((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:)=elStress((1:gl.INTPOINT),:);
    if (solverStatus=='new')
      STATEV((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:)=statev((1:gl.INTPOINT),:);
    end
    %
    %assembly
    for nodeRow=1:nElnode
      posRow = elements(e).getNodeIndex(nodeRow);
      %internal load
      Fi((posRow-1)*gl.NDOF+(1:gl.NDOF))=Fi((posRow -1)*gl.NDOF+(1:gl.NDOF))+ ...
                                       elFi((nodeRow-1)*gl.NDOF+(1:gl.NDOF));
      for nodeCol=1:nElnode
        posCol = elements(e).getNodeIndex(nodeCol);
        for i=1:gl.NDOF
          row   = gl.NDOF*(posRow -1)+i;
          rowEl = gl.NDOF*(nodeRow-1)+i;
          for m=1:gl.NDOF
            col   = gl.NDOF*(posCol -1)+m;
            colEl = gl.NDOF*(nodeCol-1)+m;
            K(row,col)=K(row,col)+elK(rowEl,colEl);
          end
        end
      end
    end
  end
end
%--------------------------------------------------------------------------
%assemblaggio dei carichi
function [F]=loadAssembly(nodes,loads,tables,Time)
  gl= GLOBAL();
  dimNode = size(nodes   ,1);
  dimLoad = size(loads   ,1);
  dimComp = gl.NDOF;
  if dimNode==0
      return;
  end
  %
  F = zeros(dimNode*gl.NDOF,1);
  %
  %assegnazione dei carichi
  for l=1:dimLoad
      indexNode = loads(l,1);
      loadInc   = Time(2);
      try
        indexTab = fix(loads(l,gl.NDOF+2));
        if indexTab~=0
          tab      = tables(indexTab);
          loadInc  = tab.getTableVal(Time(2));
        end
      catch
      end
      %
      F(gl.NDOF*(indexNode-1)+(1:dimComp))=loadInc*loads(l,(1:dimComp)+1);
  end
end
%--------------------------------------------------------------------------
function [temp]=getTemperatureVector(nodes,nodeTemp,env)
  dim = size(nodes,1);
  temp = zeros(dim,2);
  temp(1:dim)=env.TEMPERATURE;
  dim = size(nodeTemp,1);
  for i=1:dim
    index = nodeTemp(i,1);
    temp(index,1)=nodeTemp(i,2);
    try %table
      temp(index,2)=nodeTemp(i,3);
    catch
    end
  end 
end
%--------------------------------------------------------------------------
%temperature assempbly
function [T]=tempAssembly(temperatures,env,tables,Time)
  dimTemp = size(temperatures,1);
  T       = zeros(dimTemp,2);
  for i=1:dimTemp
    tempInc(1)  = Time(1);
    tempInc(2)  = Time(2);
    indexTab = temperatures(i,2);
    try %table
      tab       = tables(indexTab);
      tempInc(2)= tab.getTableVal(Time(2));
      T(i,2)=tempInc(2)*temperatures(i);
    catch
      T(i,2)=tempInc(2)*(temperatures(i)-env.TEMPERATURE)+env.TEMPERATURE;
    end
    T(i,1)=env.TEMPERATURE;
    %quando calcoliamo gli stress abbiamo bisogno della deformazione elastica
    %totale quindi qui sclcolo dT totale
  end
end
%--------------------------------------------------------------------------
function [vBC]=getBCvector(boundary,dimProblem)
  gl= GLOBAL();
  vBC=zeros(dimProblem,1);
  dimBound=size(boundary,1);
  for bc=1:dimBound
    indexNode = boundary(bc,1);
    for i=1:gl.NDOF
      if boundary(bc,i+1)==1 %il primo valore è l'indice di nodo
          indexPos      = gl.NDOF*(indexNode-1)+i;
          vBC(indexPos) = 1;
      end
    end
  end
end
%--------------------------------------------------------------------------
function [u]=solver(K,F,vBC,boundary)
  %gl= GLOBAL();
  dimProblem   = size(K,1);
  dimBound     = size(boundary,1);
  dimComp      = size(boundary,2);
  bcCont       = 0;
  indexNodeOld = 0;

  u   =zeros(dimProblem,1);

  for b=1:dimProblem
    if (vBC(b)==1)
      indexPos = b-bcCont;
      K(:,indexPos)=[];
      K(indexPos,:)=[];
      F(indexPos)  =[];
      bcCont = bcCont+1;
    end
  end

  dimk=size(K,1);
  if dimk>0
    x = K\F;
  end

  cont=1;
  for i=1:dimProblem
    if vBC(i)==0.0
      u(i) = x(cont);
      cont = cont+1;
    end
  end
end
%--------------------------------------------------------------------------
function [STRAIN]=updateSTRAIN(nodes,elements,displacements)
  gl= GLOBAL();
  dimEl = size(elements,1);
  STRAIN=zeros(dimEl*gl.INTPOINT ,     gl.SIX);
  nDisp =zeros(gl.TOTELNODES*gl.NDOF,       1);
  X     =zeros(gl.TOTELNODES     ,gl.DIMSPACE);
  for e=1:dimEl
    for n=1:gl.TOTELNODES
      nIndex     =elements(e).getNodeIndex(n);
      nDisp((n-1)*gl.NDOF+(1:gl.NDOF))=displacements((nIndex-1)*gl.NDOF+(1:gl.NDOF));
      X(n,:)                    =nodes        ( nIndex,:);
    end
    elStrain=elements(e).getElStrain(X,nDisp);
    STRAIN((e-1)*gl.INTPOINT+(1:gl.INTPOINT),:)=elStrain(1:gl.INTPOINT,:);
  end
end
%--------------------------------------------------------------------------
function [Norm]=getRedidueNorm(residue,vBC,iter)
  %gl= GLOBAL();
  Norm  = 0.0;
  r     = [];
  r     = residue;
  dimProblem = size(residue,1);
  for i=1:dimProblem
    r(i)=(1-vBC(i))*r(i);
  end
  Norm=norm(r,2);
end
%--------------------------------------------------------------------------
function [vF]=getReactions(K,u)
  %gl= GLOBAL();
  dim = size(K,1);
  vF   = zeros(dim,1);
  for i=1:dim
    vF(i)=vF(i)+K(i,1:dim)*u(1:dim);
  end
end
%--------------------------------------------------------------------------
function [var]=updateDisp(u,boundary,tables,time)
  gl=GLOBAL();
  dimBound= size(boundary,1);
  var = u;
  for bc=1:dimBound
    indexNode = boundary(bc,1);
    bcInc = time;
    try
      indexTab = int(boundary(bc,2*gl.NDOF+2));
      if indexTab~=0
        tab   = tables(indexTab);
        bcInc = tab.getTableVal(time);
      end
    catch
    end
    for dof=1:gl.NDOF
      if boundary(bc,1+dof)==1
        var(gl.NDOF*(indexNode-1)+dof)=bcInc*boundary(bc,1+gl.NDOF+dof);
      end
    end
  end
end
%--------------------------------------------------------------------------
function [var]=add_du_from_BC(u,boundary,tables,Time)
  gl=GLOBAL();
  dimBound= size(boundary,1);
  var = u;
  for bc=1:dimBound
    indexNode = boundary(bc,1);
    bcInc = Time(2)-Time(1);
    try
      indexTab = fix(boundary(bc,2*gl.NDOF+2));
      if indexTab~=0
        tab   = tables(indexTab);
        bcInc = tab.getTableVal(Time(2))-tab.getTableVal(Time(1));
      end
    catch
      %
    end
    %
    for dof=1:gl.NDOF
      if boundary(bc,1+dof)==1
        var(gl.NDOF*(indexNode-1)+dof)=bcInc*boundary(bc,1+gl.NDOF+dof);
      end
    end
  end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------