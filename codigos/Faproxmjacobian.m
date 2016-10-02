function [J]=Faproxmjacobian(Fk,p_new,p_old1,w,s,parameter,nflagface,nflagno,auxflag,gamma,weightDMP,mobility)
global elem pmethod
nelem=size(elem,1);
J=sparse(nelem,nelem);
pj=p_old1;

for ielem=1:nelem
    % atualiza a i-ésima coluna com xi+h
    pj(ielem)=p_new(ielem);
    
    % Interpolação das pressões nas faces
      % calculo da pressão
    [pinterp_new]=Fpressureinterp(pj,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
    
    if strcmp(pmethod,'nlfvDMP')
 
        % Montagem da matriz global
        
        [auxM,auxRHS]=FassemblematrixnlfvLPEW(pinterp_new,parameter,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [auxM,auxRHS] = addsource(sparse(auxM),auxRHS,wells);
        
    elseif strcmp(pmethod,'nlfvLPEW')
        
        % Montagem da matriz global
        
         [auxM,auxRHS] =FassemblematrixnlfvLPEW(pinterp_new,parameter,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [auxM,auxRHS]  = addsource(sparse(auxM),auxRHS,wells);
        
    end
    % Calculo do residou
    Fkk= auxM*pj - auxRHS;
    
    % Montagem da matriz Jacobiano por método Diferencia Finita
    J(1:nelem,ielem)=(Fkk(:)-Fk(:))./(p_new(ielem)-p_old1(ielem));
    
    % Atualiza o vetor "pj"
    pj=p_old1;
    
end

end