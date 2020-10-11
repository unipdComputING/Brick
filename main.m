%
clear;
clc;
fprintf('');
fprintf('*****************************************************\n');
fprintf('*****************************************************\n');
fprintf('***********************************************BRICK8\n');
fprintf('*****************************************************\n');
fprintf('*****************************************************\n');
fprintf('*****************************************************\n');
%inizializzazione
unitRes  =  0;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------------------------------------------------------------------INPUT
[nodes,elements,materials,boundary,loads,temperatures,...
                                tables,env,solverPar,...
                                             stFilePath]=READ_INPUT_FILE();
%
timeStart = solverPar(1);
timeEnd   = solverPar(2);
totStep   = solverPar(3);
saveEvery = solverPar(4);
%-----------------------------------------------------------------END INPUT
lenPath = length(stFilePath);
if lenPath>4 
  st        = strsplit(stFilePath,'.inp');
  stFilePath=st{1};
else
  stFilePath = 'resBRICK';
end
gid = GID();
gid.mshBEAM(stFilePath, nodes,elements);
unitRes = gid.openCloseGIDResFile(stFilePath,unitRes,'open');
%
NLINEARSOLVER(unitRes,timeStart,timeEnd,totStep,saveEvery,  ...
                               elements, nodes, materials,  ...
                               loads, temperatures,boundary,...
                               tables,env);
%
unitRes = gid.openCloseGIDResFile(stFilePath,unitRes,'close');
fprintf('END\n');
fprintf('*****************************************************\n');
fprintf('*****************************************************\n');
fprintf('*****************************************************\n');



