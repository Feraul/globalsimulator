function  [p,flowrate,flowresult]=FsolverPressurelfvHP(mobility,wells,parameter,weightDMP,nflagface,nflagno,auxflag)
%Initialize "bedgesize" and "inedgesize"
w=0;
s=0;
% Montagem da matriz global
 [M,I]=FassemblematrixlfvHP(parameter,nflagface,weightDMP,mobility);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p =M\I;

%Message to user:
disp('>> The Pressure field was calculated with success!');

[pinterp]=Fpressureinterp(p,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);

%Get the flow rate 
[flowrate,flowresult]=FflowratelfvHP(parameter,weightDMP,mobility,pinterp,p);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
