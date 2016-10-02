function [p,flowrate,flowresult]=FsolverPressurelfvLPEW(kmap,nflagface,...
    nflagno,mobility,wells,Sw,parameter,weightDMP,V,N,auxflag)

global interptype phasekey bedge inedge 


%It switches according to "interptype"
switch char(interptype)
    %LPEW 1
    case 'lpew1'
        % calculo dos pesos que correspondem ao LPEW1
        [w,s] = ferncodes_Pre_LPEW_1(kmap,mobility,V,Sw,N);
    %LPEW 2
    case 'lpew2'
        % calculo dos pesos que correspondem ao LPEW2
        [w,s] = ferncodes_Pre_LPEW_2(kmap,mobility,V,Sw,N);
end  %End of SWITCH

% Montagem da matriz global
[M,I]=FassemblematrixlfvLPEW(parameter,w,s,nflagno,weightDMP,mobility);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p = M\I;

%Message to user:
disp('>> The Pressure field was calculated with success!');

[pinterp]=Fpressureinterp(p,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
%Get the flow rate (Diamond)
[flowrate,flowresult]=FflowratelfvLPEW(parameter,weightDMP,mobility,pinterp,p);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
