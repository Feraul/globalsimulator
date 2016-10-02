
function [p,step,errorelativo,flowrate,flowresult] = JFNK1(tolnewton,kmap,...
                parameter,metodoP,auxflag,w,s,nflag,fonte,gamma,...
                nflagno,benchmark,M_old1,RHS_old1,p_old1,R0,weightDMP,auxface)
% inicializando
R= M_old1*p_old1-RHS_old1;
er=1;
step= 0; % iteration counter
while tolnewton<er
    step = step + 1; % update iteration counter
   
    % O método de Newton- Krylov aqui implementado foi proposto no
    % artigo "Jacobian-free Newton–Krylov methods:a survey of approaches and applications
    % D.A. Knoll, D.E. Keyes e foi implementado no Matlab segui o site:
    % http://www.mathworks.com/matlabcentral/fileexchange/45170-jacobian-free-newton-krylov--jfnk--method
    % para sistema de equações dadas, nós adaptamos a nosso contexto.
    j_v_approx1 = @(v)JV_APPROX1(v, R,p_old1,nflag,w,s,metodoP,parameter,...
        kmap,nflagno,benchmark,fonte,auxflag,gamma,weightDMP,auxface);
    
    % Calculo pelo método de iteração GMRES do Matlab
    %[v,flag,relres,iter,resvec] = gmres(j_v_approx1, R,2,1e-5,size(elem,1)*0.5); % solve for Krylov vector
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    % Calculo pelo método de iteração BICGSTAB do Matlab
    [v,flag]=bicgstabl(j_v_approx1,R,1e-6,200);
    flag
    
    % Calculo da pressão na iteração k+1
    p_new = p_old1 - v;
    % garantido a positividade
    n=1;
    while min(p_new)<0
        p_new=p_old1-v*(1/2^n);
        n=2*n;
    end
    
    % Plotagem no visit
    postprocessor(p_new,step)
    
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,w,s,auxflag,metodoP,parameter,weightDMP);
    
    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflag,nflagno...
                ,parameter,kmap,fonte,metodoP,w,benchmark,weightDMP,auxface);
    
    % Calculo do residuo
    R= M_new*p_new - RHS_new;
    
    % calculo do erro
    A=logical(norm(R0) ~= 0.0);
    B=logical(norm(R0) == 0.0);
    er=A*abs(norm(R)/norm(R0))+B*0
    errorelativo(step)=er;
    
    % Atualizando a pressão
    p_old1=p_new;
    
end
p=p_old1;
pinterp=pressureinterp(p,nflag,w,s,auxflag,metodoP,parameter,weightDMP);
if strcmp(metodoP,'nlfvDMPSY')
% implementação do fluxo NLFV-DMP
[flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflag,kmap,gamma,weightDMP);
else
% implementação do fluxo NLFV
[flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter);
end
end
