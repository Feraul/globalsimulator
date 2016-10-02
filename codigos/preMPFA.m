%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 31/07/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the 

%--------------------------------------------------------------------------

function [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
    Bleft,Bright,overedgecoord,normk,Fg,mapinv,maptransm,mapknownvec,...
    pointedge,bodyterm,Hesq,Kde,Kn,Kt,Ded,V,N,kmap,nflag,parameter,gamma,...
    tolpicard,niter,facelement,nflagface,weightDMP] = preMPFA(kmap,...
    klb)
%Define global parameters:
global pmethod;

%Obtain the coordinate of both CENTER and AUXILARY nodes of elements which
%constitute the mash. The AREA of each element is also calculated.
%Message to user:
disp(' ');
disp('---------------------------------------------------');
disp('>> Preprocessing Pressure Equation...');
disp(' ');

%Initialize all parameters of output
transmvecleft = 0;
transmvecright = 0;
knownvecleft = 0;
knownvecright = 0;
storeinv = 0;
Bleft = 0;
Bright = 0;
Fg = 0;
mapinv = 0;
maptransm = 0;
mapknownvec = 0;
pointedge = 0;
bodyterm = 0;
Hesq = 0;
Kde = 0;
Kn = 0;
Kt = 0;
Ded = 0;
V = 0;
N = 0;
nflag = 0;
nflagface=0;
weightDMP=0;
parameter=0;
gamma=0;
tolpicard=0;
niter=0;
facelement=0;
%Define parametric variables:
%Parameter Used in Full Pressure Support (FPS) 
%"p" quadrature point to flux in the auxilary sub interaction region
p = 1;
%Parameter Used in Full Pressure Support (FPS) and Triangle Pressure 
%Support (TPS)  
%"q" quadrature point to flux in the sub interaction region
q = 1;

%Fill the matrix "overedgecoord"
overedgecoord = overedge;
%Define the norm of permeability tensor ("normk")
[normk,kmap] = calcnormk(kmap);

%Get the length of the edge with non-null Neumann Boundary Condition.
knownboundlength = getknownboundlength(klb);

%--------------------------------------------------------------------------
%Calculate the TRANSMISSIBILITY parameters:

%Chose the type of MPFA according "pmethod"
switch char(pmethod)
    %Calculate the transmissibilities from TPFA
    case 'tpfa'
        [transmvecleft,knownvecleft,Fg,bodyterm] = transmTPFA(kmap);    
    %Calculate the little matrices for MPFA-TPS (Aavatsmark et al., 1998)
    case 'mpfao'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmTPS(kmap,q,knownboundlength);
    %Calculate the little matrices for MPFA-FPS (Edwards and Zheng, 2008)
    case 'fps'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmFPS(kmap,q,p,knownboundlength);
    %Calculate the little matrices for MPFA-Enriched (Chen et al., 2008)
    case 'empfa'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmEnriched(kmap,knownboundlength);    
    %Calculate geometrical and physical terms to be used in MPFA-Diamond 
    %(Gao and Wu, 2010). 
    case 'mpfad'
        %Get "ferncodes_nflag" 
        nflag = ferncodes_calflag;
        %Get preprocessed terms:
        [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap);
        %Call another parameters that I don't know.
        [V,N,] = ferncodes_elementface; 
    case 'mpfahd'
        %% calculo dos parametros ou constantes (ksi)
        [parameter]=FM2coefficient(kmap); % baseado nos senos
        %[parameter]=FM3coefficient(kmap); % baseado na teoria do co-normal
        %[parameter]=FcoefficientLPSangle(kmap);
        % adequação dos flags de contorno
        [nflag]= Fcontflagno;
        % calculo dos pesos DMP
        [weightDMP]= FweightnlfvDMP(kmap);
        %Call another parameters that I don't know.
        [V,N,~] = ferncodes_elementface;
        
    case {'mpfah','nlfvdmp'}
        % faces alrededor de um elemento
        [V,N,facelement]=ferncodes_elementface;
        % calculoa dos pontos armonicos
        [pointarmonic]=Fharmonicopoint(kmap,N);
        % calculo dos parametros ou constantes (ksi)
        %[parameter]=FcoefficientPPSharmonicpoint(facelement,pointarmonic,kmap);
        [parameter]=FM1coefficient(facelement,pointarmonic,kmap); % baseado no produto misto
        % adequação dos flag de face de contorno
        [nflagface]= Fcontflagface;
        % adequação dos nos flags de contorno
        [nflag]= Fcontflagno;
        % calculo dos pesos DMP
        [weightDMP]= FweightnlfvDMP(kmap);    
        % dados inicialização métodos dos volumes finitos não linear
        gamma=0.0; % este parametro esta no intervalo [0,1]
        tolpicard=1e-07; % tolerancia
        niter=1000; % iterações
        
    case 'nlfvpp'
        % calculo dos parametros ou constantes (ksi)
        
        
        %[parameter]=FcoefficientLPSangle(kmap); %trabajado com esto no primeiro teste
        [parameter]=FM2coefficient(kmap); % baseado nos senos
        %[parameter]=FM3coefficient(kmap); % baseado na teoria do co-normal
        % adequação dos flag de face de contorno
        [nflagface]= Fcontflagface;
        % adequação dos nos flags de contorno
        [nflag]= Fcontflagno;
        % calculo dos pesos DMP
        [weightDMP]= FweightnlfvDMP(kmap);
        %Call another parameters that I don't know.
        [V,N,~] = ferncodes_elementface;
        % dados inicialização métodos dos volumes finitos não linear
        gamma=0.0; % este parametro esta no intervalo [0,1]
        tolpicard=1e-07; % tolerancia
        niter=1000; % iterações
end  %End of SWITCH

%Message to user:
disp('>> "preMPFA" was finished with success!');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "overedgesection"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordsection] = overedgesection(edgematrix)
%Define global parameters:
global coord;

%Initialize the matrix used in this function
coordsection = zeros(size(edgematrix,1),3);
%This loop swept among the limits stabelished above
%"firstlim" and "lastlim" are parameters which define where the loop begins
%and where its finish

iover = 1:size(edgematrix,1);
coordsection(iover,1:3) = 0.5*(coord(edgematrix(iover,1),:) + ...
    coord(edgematrix(iover,2),:));

%--------------------------------------------------------------------------
%Function "overedge"
%--------------------------------------------------------------------------

function [overedgecoord] = overedge
%Define global parameters:
global bedge inedge;

%Initialize the matrix "overedgecoord"
overedgecoord = zeros((size(bedge,1) + size(inedge,1)),3);

%Fill "overedgecoord" (just edge over boundary)
overedgecoord(1:size(bedge,1),:) = overedgesection(bedge);
%Fill the "overedgecoord" rest (just edge inside domain)
%"continedge" is an internal edge's counter
overedgecoord(size(bedge,1) + 1:size(overedgecoord,1),:) = ...
    overedgesection(inedge);

%--------------------------------------------------------------------------
%Function "calcnormk"
%--------------------------------------------------------------------------

function [normk,kmap] = calcnormk(kmap)
%Define global parameters:
global elem centelem;

%Initialize "normk" (it is a vector)
normk = zeros(size(centelem,1),1);
%Define the norm of permeability tensor
%Obtain "kmap" for each case
kmap = PLUG_kfunction(kmap);
%Swept all elements
for ik = 1:length(normk)
    %Define the material pointer in "elem"
    pointer = elem(ik,5);
    %It catches only the permeability components
    permcompon = [kmap(pointer,2) kmap(pointer,3); ...
        kmap(pointer,4) kmap(pointer,5)];
    %Calculate the norm of tensor
    normk(ik) = norm(permcompon);
end  %End of FOR
