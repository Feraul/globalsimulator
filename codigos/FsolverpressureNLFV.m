function [p,flowrate,flowresult]=FsolverpressureNLFV(kmap,nflagface,...
    nflagno,mobility,Sw,parameter,auxflag,niter,tol,iteration,weightDMP,gamma,V,N,wells)
%Define global parameters
global interptype elem pmethod 

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
% pressure inicialization
p_old=0.1*ones(size(elem,1),1);
% interpolação nos nós ou faces
[pinterp]=Fpressureinterp(p_old,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
% Montagem da matriz global
[M,I]=FmatrixNLFV(p_old,pinterp,parameter,weightDMP,mobility,wells);

if strcmp(iteration,'iterpicard')
    
    [p]=Fiterpicard(M,I,niter,tol,parameter,auxflag,w,s,nflagface,nflagno,p_old,weightDMP,mobility,wells);
    
elseif strcmp(iteration, 'iterdiscretnewton')
    
    p_old1=M\I;
    % interpolação nos nós ou faces
    [pinterp1]=Fpressureinterp(p_old,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
    
    % montagem global matriz
    [M1,I1]=FmatrixNLFV(pinterp1,parameter,mobility,wells);
    
    % resolvedor de pressão pelo método de Newton-Discreto
    [p]=Fiterdiscretnewton(M1,I1,tol,nflagface,nflagno,...
        parameter,auxflag,w,s,p_old,gamma,M,I,p_old1,weightDMP,mobility,wells);
    
    
elseif strcmp(iteration, 'JFNK')
    
    
    p_old1=M\I;
    % calculo do residuo
    R0=M*p_old-I;
    
    % interpolação nos nós ou faces
    [pinterp1]=Fpressureinterp(p_old1,nflagface,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
    
    % montagem global matriz
    [M1,I1]=FmatrixNLFV(pinterp1,parameter,mobility,wells);
    
    % calculo da pressão
    [p]= FJFNK1(tol,kmap,parameter,metodoP,auxflag,w,s,nflagface,fonte,gamma,nflagno,M1,I1,p_old1,R0,weightDMP,auxface,mobility);
    
end

%Message to user:
disp('>> The Pressure field was calculated with success!');

[pinterp]=Fpressureinterp(p,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);

if strcmp(pmethod,'nlfvpp')
    %Get the flow rate
    [flowrate,flowresult]=FflowrateNLFV(p, pinterp, parameter,mobility);
    
else
 % 
 [flowrate,flowresult]=FflowrateNLFVDMP(p, pinterp, parameter,nflagface,kmap,gamma,weightDMP,mobility);
end

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
