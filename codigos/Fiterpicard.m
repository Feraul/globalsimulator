function [p]=Fiterpicard(M_old,RHS_old,niter,tol,parameter,auxflag,...
                         w,s,nflagface,nflagno,p_old,weightDMP,mobility,wells)

% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

% inicializando dados para iteração Picard
step=0;
er=1;
while (tol<er || tol==er) && (step<niter)
    % atualiza iterações
    step=step+1;
    %% calculo da pressão utilizando o precondicionador
    
    %[L,U] = ilu(M_old,struct('type','crout','droptol',1e-6));
    %M=L*U;
    
    %u=(M_old*inv(M))\RHS_old;
    
    %p_new=M\u;
    % calculo da pressão pelo método direto
    %p_new=inv(M_old)*RHS_old;  % inversão com pivotamento
    
    p_new=M_old\RHS_old;  % inversão sem pivotamento
    
    % calculo da pressão usando método iterativo + precondicionador ILU
    %[L,U]=ilu(M_old,struct('type','ilutp','droptol',1e-6));
    %M=L*U;
    % GMRES do artigo SIAM
    %[p_new, error, iter, flag] = gmresSIAM( M_old, p_old, RHS_old, M, 10, 1000, 1e-10 );
    
    % GMRES do MATLAB
    %[p_new,flag]=gmres(M_old,RHS_old,2,1e-6,100,L,U); % gmres Matlab
    %flag
    
    % plotagem no visit
    Fpostprocessor(p_new,step)
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=Fpressureinterp(p_new,nflagface,nflagno,w,s,auxflag,parameter,weightDMP,mobility);
   
     % montagem global matriz
    [M_new,RHS_new]=FmatrixNLFV(p_new,pinterp_new,parameter,weightDMP,mobility,wells);
    
    %% Calculo do residuo
    R = norm(M_new*p_new - RHS_new);
    
    %% calculo do erro
    A=logical(R0 ~= 0.0);
    B=logical(R0 == 0.0);
    er=A*abs(R/R0)+B*0;
    
    %% atualizar
    M_old=M_new;
    RHS_old=RHS_new;
    
end
p=M_old\RHS_old;

residuo=er;
niteracoes=step;

disp('>> The Pressure field is calculated with "nlfvpp"')

x=['>> Tolerance reached:',num2str(residuo)];
disp(x);
y=['>> Number of iterations:',num2str(niteracoes)];
disp(y);
end