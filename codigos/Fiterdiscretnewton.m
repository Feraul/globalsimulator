function [p]=Fiterdiscretnewton(M_old,RHS_old,tol,nflagface,nflagno,...
           parameter,auxflag,w,s,p_old,gamma,M_old1,RHS_old1,p_old1,weightDMP,mobility)

global pmethod

% inicializando dados para iteração Picard
Rk= M_old1*p_old1-RHS_old1;
R0= M_old*p_old-RHS_old;
er=1;
step=0;

% calculo do jacobiano discreto
[J]=Faproxmjacobian(R0,p_old1,p_oldw,s,parameter,nflagface,nflagno,auxflag,gamma,weightDMP,mobility);
%
while tol<er
    step=step+1;
    
    % calculo inversa da Jaconiano
    pr=-J\Rk;
    
    % calculo da pressão
    p_new=p_old1+pr;
    
    % plotagem no visit
    %postprocessor(p_new,step)
    
    % calculo da pressão
    [pinterp_new]=Fpressureinterp(p_new,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
    
    if strcmp(pmethod,'nlfvDMP')
 
        % Montagem da matriz global
        
        [M_new,RHS_new]=FassemblematrixnlfvLPEW(pinterp_new,parameter,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [M_new,RHS_new] = addsource(sparse(M_new),RHS_new,wells);
        
    elseif strcmp(pmethod,'nlfvLPEW')
        
        % Montagem da matriz global
        
         [M_new,RHS_new] =FassemblematrixnlfvLPEW(pinterp_new,parameter,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [M_new,RHS_new]  = addsource(sparse(M_new),RHS_new,wells);
        
    end
    
    % calculo do residuo
    Rkk= M_new*p_new - RHS_new;
    
    % calculo do Jacobiano discreto
    [J]=Faproxmjacobian(Rk,p_new,p_old1,w,s,parameter,nflagface,nflagno,auxflag,gamma,weightDMP,mobility);
    
    % calculo do erro
    A=logical(norm(R0) ~= 0.0);
    B=logical(norm(R0) == 0.0);
    er=A*abs(norm(Rkk)/norm(R0))+B*0
    errorelativo(step)=er;
    
    % atualizar
    p_old1=p_new;
    Rk=Rkk;
end
p=p_old1;
end