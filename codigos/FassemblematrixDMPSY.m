function [M,I]=FassemblematrixDMPSY(p,pinterp,gamma,parameter,weightDMP,mobility)
global inedge coord bedge bcflag elem phasekey
bedgesize = size(bedge,1);
valuemin=1e-16;
I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));

for ifacont=1:size(bedge,1)
    
    %Define "mobonface" (for "bedge")
    %It is a One-phase flow. In this case, "mobility" is "1"
    if phasekey == 1
        mobonface = mobility;
    %It is a Two-phase flow
    else
        %"mobonface" receivees the summation of water and oil
        %mobilities (came from "IMPES" - Marcio's code modification)
        mobonface = sum(mobility(ifacont,:));
    end  %End of IF
    % atirbui o elemento a esquerda da face IJ
    lef=bedge(ifacont,3);
    % calcula a norma da face IJ
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        % atribui os flag do fluxo prescrito 
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        % atribui o valor do fluxo prescrito
        I(lef)=I(lef)- normcont*bcflag(r,2);
    else
        % calculo da mobilidade
        auxmobility=mobonface;
        % faces que conformam os eixos auxiliares
        ifacelef1= parameter(1,3,ifacont);
        % montagem da matriz global correspondente à face "ifacelef1"
        [M,I]=FauxassemblematrixcontourDMPSY(ifacelef1,M,I,parameter(1,1,ifacont),lef,pinterp,normcont,weightDMP,auxmobility);
        % faces que conformam os eixos auxiliares
        ifacelef2= parameter(1,4,ifacont);
        % montagem da matriz global correspondente à face "ifacelef2"
        [M,I]=FauxassemblematrixcontourDMPSY(ifacelef2,M,I,parameter(1,2,ifacont),lef,pinterp,normcont,weightDMP,auxmobility);
        
    end
    
end

%% Montagem da matriz global
x=1000;
y=1000;
for iface=1:size(inedge,1)
    
    
    %Define "mobonface" (for "bedge")
    %It is a One-phase flow. In this case, "mobility" is "1"
    if phasekey == 1
        mobonface = mobility;
    %It is a Two-phase flow
    else
        %"mobonface" receivees the summation of water and oil
        %mobilities (came from "IMPES" - Marcio's code modification)
        mobonface = sum(mobility(bedgesize + iface,:));
    end  %End of IF
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    
    %% Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    % calculo da mobilidade
    auxmobility=mobonface;
    % Fluxo parcial do elemento a esquerda
    
    F1=auxmobility*norma*(parameter(1,1,ifactual)*(p(lef)-pinterp(parameter(1,3,ifactual)))+ parameter(1,2,ifactual)*(p(lef)-pinterp(parameter(1,4,ifactual))));
    
    % Fluxo parcial do elemento a direita
    F2=auxmobility*norma*(parameter(2,1,ifactual)*(p(rel)-pinterp(parameter(2,3,ifactual)))+ parameter(2,2,ifactual)*(p(rel)-pinterp(parameter(2,4,ifactual))));
    if abs(F1)<1e-20
        F1=0;
    end
    if abs(F2)<1e-20
        F2=0;
    end
    %% quando F1*F2<=0
    if F2*F1<0 || F2*F1==0
        % Artigo de Gao and Wu, 2013
        
        mu1=(abs(F2)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        mu2=(abs(F1)+valuemin)/(abs(F1)+abs(F2)+2*valuemin);
        
        betalef=mu1*(1-sign(F1*F2));
        betarel=mu2*(1-sign(F1*F2));
        
        % contribuições do elemento a esquerda
        
        ifacelef1=parameter(1,3,ifactual);
        [M,I]=FauxassemblematrixintDMPSY1(ifacelef1,M,I,parameter(1,1,ifactual),lef,pinterp,norma,betalef,weightDMP,auxmobility);
        
        ifacelef2=parameter(1,4,ifactual);
        [M,I]=FauxassemblematrixintDMPSY1(ifacelef2,M,I,parameter(1,2,ifactual),lef,pinterp,norma,betalef,weightDMP,auxmobility);
        
        % contribuições do elemento a direita
        
        ifacerel1= parameter(2,3,ifactual);
        [M,I]=FauxassemblematrixintDMPSY1(ifacerel1,M,I,parameter(2,1,ifactual),rel,pinterp,norma,betarel,weightDMP,auxmobility);
        
        ifacerel2= parameter(2,4,ifactual);
        [M,I]=FauxassemblematrixintDMPSY1(ifacerel2,M,I,parameter(2,2,ifactual),rel,pinterp,norma,betarel,weightDMP,auxmobility);
        
    else %% F1*F2>0
        [F1b,]=Fcalfluxopartial2(parameter(1,3,ifactual),parameter(1,4,ifactual),...
            parameter(1,1,ifactual), parameter(1,2,ifactual), gamma,lef,pinterp,ifactual,p,norma,auxmobility);
        % aveliando no elemento a direita
        
        [F2b,]=Fcalfluxopartial2(parameter(2,3,ifactual),parameter(2,4,ifactual),...
            parameter(2,1,ifactual), parameter(2,2,ifactual), gamma,rel,pinterp,ifactual,p,norma,auxmobility);
        
        %% Calculo das constantes da não linearidade (eq. 13) do artigo Gao e Wu, (2013)
        % veja a REMARK 3.2 pag. 313 do mesmo artigo
        if F1b*F2b>0
            mulef=abs(F2b)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=abs(F1b)/(abs(F1b)+abs(F2b)+2*valuemin);
        else
            mulef=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
            
            murel=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        end
        
        
        mulef1=(abs(F2b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        
        murel1=(abs(F1b)+valuemin)/(abs(F1b)+abs(F2b)+2*valuemin);
        %
        ifacelef1=parameter(1,3,ifactual);
        ifacelef2=parameter(1,4,ifactual);
        ifacerel1= parameter(2,3,ifactual);
        ifacerel2= parameter(2,4,ifactual);
        if ifacerel1==ifactual
            auxparameter2=parameter(2,1,ifactual);
        elseif ifacerel2==ifactual
            
            auxparameter2=parameter(2,2,ifactual);
        else
            auxparameter2=0;
            x=0;
            
        end
        %------------------------------------------------------------%
        if ifacelef1==ifactual
            auxparameter1=parameter(1,1,ifactual);
        elseif ifacelef2==ifactual
            
            auxparameter1=parameter(1,2,ifactual);
        else
            auxparameter1=0;
            y=0;
            
        end
        %% implementação do beta
        betalef=mulef*(1-sign(F1b*F2b));
        betarel=murel*(1-sign(F1b*F2b));
        %=================================================================%
        if F1b*F2b<0 || F1b*F2b==0 || (F2b<0 && F1b<0 && x==0 && y==0) || (F2b>0 && F1b>0 && x==0 && y==0)
            %% Calculo das contribuições do elemento a esquerda
            % [weightlef,weightrel]=weightnlfvDMP(kmap,ifactual-size(bedge,1));
            % Calculo das contribuições do elemento a esquerda
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            
            %% contribuição da transmisibilidade no elemento esquerda
            M(lef,lef)=M(lef,lef)+ alfa*auxmobility;
            M(lef,rel)=M(lef,rel)- alfa*auxmobility;
            
            [M,I]=Fauxassemblematrixint(ifacelef1,M,I,parameter(1,1,ifactual),auxparameter2,...
                lef,rel,mulef1,murel1,weightlef,weightrel,pinterp,norma,betalef,ifactual,gamma,weightDMP,auxmobility);
            
            [M,I]=Fauxassemblematrixint(ifacelef2,M,I,parameter(1,2,ifactual),auxparameter2,...
                lef,rel,mulef1,murel1,weightlef,weightrel,pinterp,norma,betalef,ifactual,gamma,weightDMP,auxmobility);
            
            %% contribuição da transmisibilidade no elemento direita
            
            M(rel,rel)=M(rel,rel)+ alfa*auxmobility;
            M(rel,lef)=M(rel,lef)- alfa*auxmobility;
            
            [M,I]=Fauxassemblematrixint(ifacerel1,M,I,parameter(2,1,ifactual),auxparameter1,...
                rel,lef,murel1,mulef1,weightrel,weightlef,pinterp,norma,betarel,ifactual,gamma,weightDMP,auxmobility);
            
            
            [M,I]=Fauxassemblematrixint(ifacerel2,M,I,parameter(2,2,ifactual),auxparameter1,...
                rel,lef,murel1,mulef1,weightrel,weightlef,pinterp,norma,betarel,ifactual,gamma,weightDMP,auxmobility);
            %=============================================================%
            %% F1b*F2b>0
        else
            
            %[weightlef,weightrel]=weightnlfvDMP(kmap,ifactual-size(bedge,1));
            % Calculo das contribuições do elemento a esquerda
            weightlef=weightDMP(ifactual-size(bedge,1),1);
            weightrel=weightDMP(ifactual-size(bedge,1),2);
            
            % Calculo da transmisibilidade (eq. 18)
            alfa=(1-gamma)*norma*(mulef1*auxparameter1*weightrel + murel1*auxparameter2*weightlef);
            % contribuição da transmisibilidade no elemento esquerda
            M(lef,lef)=M(lef,lef)+ alfa*auxmobility;
            M(lef,rel)=M(lef,rel)- alfa*auxmobility;
            % contribuição da transmisibilidade no elemento direita
            M(rel,rel)=M(rel,rel)+ alfa*auxmobility;
            M(rel,lef)=M(rel,lef)- alfa*auxmobility;
        end
        
    end
    
end
% switch benchmark
%     case 'gaowu6'
%         %% malha 23x23
%         % M(357,:)=0*M(357,:);
%         % M(357,357)=1;
%         % I(357)=1;
%         % M(173,:)=0*M(173,:);
%         % M(173,173)=1;
%         % I(173)=0;
%         %% malha 11x11
%         M(83,:)=0*M(83,:);
%         M(83,83)=1;
%         I(83)=1;
%         M(39,:)=0*M(39,:);
%         M(39,39)=1;
%         I(39)=0;
% end
end