function [M,I]=FassemblematrixlfvLPEW(parameter,w,s,nflag,weightDMP,mobility)
global inedge coord bedge bcflag elem esurn1 esurn2 phasekey numcase;
 if numcase==31.2
        auxbcflag=bcflag;
        auxbcflag(2,2)=-1;
    else
        auxbcflag=bcflag;
 end
  
% incialização das matrizes
I=zeros(size(elem,1),1);
M=zeros(size(elem,1),size(elem,1));

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);

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
    
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)- normcont*auxbcflag(r,2);
    else
        % implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ mobonface*normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        nolef1=parameter(1,3,ifacont);
        nolef2=parameter(1,4,ifacont);
        % contribuições do nó 1
        if nflag(nolef1,1)<200
            I(lef,1)=I(lef,1)+ mobonface*normcont*(parameter(1,1,ifacont)*nflag(nolef1,2));
        else
            auxvetor=1:(esurn2(nolef1+1)-esurn2(nolef1));
            auxvetor=auxvetor+esurn2(nolef1);
            auxesurn1=esurn1(auxvetor);
            
            M(lef, auxesurn1)=M(lef, auxesurn1)-mobonface*normcont*parameter(1,1,ifacont)*w(auxvetor);
            
        end
        % contribuições do nó 2
        if nflag(nolef2,1)<200
            
            I(lef,1)=I(lef,1)+ mobonface*normcont*parameter(1,2,ifacont)*nflag(nolef2,2);
            
        else
            auxvetor=1:(esurn2(nolef2+1)-esurn2(nolef2));
            auxvetor=auxvetor+esurn2(nolef2);
            auxesurn1=esurn1(auxvetor);
            M(lef, auxesurn1)=M(lef, auxesurn1)-mobonface*normcont*parameter(1,2,ifacont)*w(auxvetor);
            
        end
        
    end
end

% Montagem da matriz global

for iface=1:size(inedge,1)
    %Define "mobonface" (for "inedge")
    %It is a One-phase flow
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
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    %[mulef, murel]=weightnlfvDMP(kmap,iface);
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxnolef1=parameter(1,3,ifactual);
    auxnolef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxnorel1=parameter(2,3,ifactual);
    auxnorel2=parameter(2,4,ifactual);
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*murel*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ARR=norma*mulef*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ mobonface*ALL;
    M(lef,rel)=M(lef,rel)- mobonface*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ mobonface*ARR;
    M(rel,lef)=M(rel,lef)- mobonface*ALL;
    % contribuições esquerda
    
    if nflag(auxnolef1,1)>200
        auxvetor=1:(esurn2(auxnolef1+1)-esurn2(auxnolef1));
        auxvetor=auxvetor+esurn2(auxnolef1);
        auxesurn1=esurn1(auxvetor);
        
        M(lef, auxesurn1)=M(lef, auxesurn1)-mobonface*murel*norma*parameter(1,1,ifactual)*w(auxvetor);
        M(rel, auxesurn1)=M(rel, auxesurn1)+mobonface*murel*norma*parameter(1,1,ifactual)*w(auxvetor);
        
    else
        I(lef,1)=I(lef,1)+ mobonface*murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
        I(rel,1)=I(rel,1)- mobonface*murel*norma*parameter(1,1,ifactual)*nflag(auxnolef1,2);
    end
    
    if nflag(auxnolef1,1)==202
        
        I(lef)=I(lef)+mobonface*murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
        
        I(rel)=I(rel)-mobonface*murel*norma*parameter(1,1,ifactual)*s(auxnolef1); %ok
    end
    
    if nflag(auxnolef2,1)>200
        
        auxvetor=1:(esurn2(auxnolef2+1)-esurn2(auxnolef2));
        auxvetor=auxvetor+esurn2(auxnolef2);
        auxesurn1=esurn1(auxvetor);
        
        M(lef, auxesurn1)=M(lef, auxesurn1)- mobonface*murel*norma*parameter(1,2,ifactual)*w(auxvetor);
        M(rel, auxesurn1)=M(rel, auxesurn1)+ mobonface*murel*norma*parameter(1,2,ifactual)*w(auxvetor);
        
    else
        I(lef,1)=I(lef,1)+ mobonface*murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        I(rel,1)=I(rel,1)- mobonface*murel*norma*parameter(1,2,ifactual)*nflag(auxnolef2,2);
        
    end
    if nflag(auxnolef2,1)==202
        
        I(lef)=I(lef)+mobonface*murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
        
        I(rel)=I(rel)-mobonface*murel*norma*parameter(1,2,ifactual)*s(auxnolef2); %ok
    end
    % direita
    if nflag(auxnorel1,1)>200
        
        
        auxvetor=1:(esurn2(auxnorel1+1)-esurn2(auxnorel1));
        auxvetor=auxvetor+esurn2(auxnorel1);
        auxesurn1=esurn1(auxvetor);
        
        M(lef, auxesurn1)=M(lef, auxesurn1)+mobonface*mulef*norma*parameter(2,1,ifactual)*w(auxvetor);
        M(rel, auxesurn1)=M(rel, auxesurn1)-mobonface*mulef*norma*parameter(2,1,ifactual)*w(auxvetor);
        
    else
        I(lef,1)=I(lef,1)- mobonface*mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
        I(rel,1)=I(rel,1)+ mobonface*mulef*norma*parameter(2,1,ifactual)*nflag(auxnorel1,2);
    end
    if nflag(auxnorel1,1)==202
        
        I(lef)=I(lef)-mobonface*mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
        
        I(rel)=I(rel)+mobonface*mulef*norma*parameter(2,1,ifactual)*s(auxnorel1); %ok
    end
    if nflag(auxnorel2,1)>200
        
        auxvetor=1:(esurn2(auxnorel2+1)-esurn2(auxnorel2));
        auxvetor=auxvetor+esurn2(auxnorel2);
        auxesurn1=esurn1(auxvetor);
        
        M(lef, auxesurn1)=M(lef, auxesurn1)+ mobonface*mulef*norma*parameter(2,2,ifactual)*w(auxvetor);
        M(rel, auxesurn1)=M(rel, auxesurn1)- mobonface*mulef*norma*parameter(2,2,ifactual)*w(auxvetor);
        
    else
        I(lef,1)=I(lef,1)- mobonface*mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
        I(rel,1)=I(rel,1)+ mobonface*mulef*norma*parameter(2,2,ifactual)*nflag(auxnorel2,2);
    end
    if nflag(auxnorel2,1)==202
        
        I(lef)=I(lef)-mobonface*mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
        
        I(rel)=I(rel)+mobonface*mulef*norma*parameter(2,2,ifactual)*s(auxnorel2); %ok
    end
end
%% malha 23x23
% M(357,:)=0*M(357,:);
% M(357,357)=1;
% I(357)=1;
% M(173,:)=0*M(173,:);
% M(173,173)=1;
% I(173)=0;
%% malha 11x11
% M(83,:)=0*M(83,:);
% M(83,83)=1;
% I(83)=1;
% M(39,:)=0*M(39,:);
% M(39,39)=1;
% I(39)=0;

end