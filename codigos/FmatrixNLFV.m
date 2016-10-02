function [M,I]=FmatrixNLFV(p,pinterp,parameter,weightDMP,mobility,wells)

global pmethod 

if strcmp(pmethod,'nlfvdmp')
    
    % Montagem da matriz global
    
    [M,I]=FassemblematrixDMPSY(p,pinterp,0,parameter,weightDMP,mobility);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    
    %Often it may change the global matrix "M"
    [M,I] = addsource(sparse(M),I,wells);
elseif strcmp(pmethod,'nlfvpp')
    
    % Montagem da matriz global
    
    [M,I] =FassemblematrixnlfvLPEW(pinterp,parameter,mobility);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    
    %Often it may change the global matrix "M"
    [M,I]  = addsource(sparse(M),I,wells);
    
end

end