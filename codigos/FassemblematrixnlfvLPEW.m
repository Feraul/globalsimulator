function [M,I]=FassemblematrixnlfvLPEW(pinterp,parameter,mobility)
global inedge coord bedge bcflag elem phasekey numcase
if numcase==31.2
        auxbcflag=bcflag;
        auxbcflag(2,2)=-1;
    else
        auxbcflag=bcflag;
 end
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
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
    
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(lef)=I(lef)- normcont*auxbcflag(r,2);
    else
        %% calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        
        alef= mobility(ifacont)*normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
            parameter(1,2,ifacont)*pinterp(parameter(1,4,ifacont)));
        
        Alef=mobility(ifacont)*normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        
        %% implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ Alef;
        I(lef,1)=I(lef,1)+alef;
    end
end

%% Montagem da matriz global

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
    norma= sqrt(vd1(1,1)^2+vd1(2)^2);
    ifactual=iface+size(bedge,1);
    
    % calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    alef=norma*(parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual)));
    
    % direita
    arel= norma*(parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual)));
    
    mulef=(abs(arel)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
    
    murel=(abs(alef)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
    
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    
    ARR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ mobonface*ALL;
    M(lef,rel)=M(lef,rel)- mobonface*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ mobonface*ARR;
    M(rel,lef)=M(rel,lef)- mobonface*ALL;
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