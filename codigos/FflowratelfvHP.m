function [flowrate,flowresult]=FflowratelfvHP(parameter,weightDMP,mobility,pinterp,p)
global inedge coord bedge bcflag centelem phasekey

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
    
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
        flowrate(ifacont,1)=normcont*bcflag(r,2);
    else
        facelef1=parameter(1,3,ifacont);
        facelef2=parameter(1,4,ifacont);
        
        flowrate(ifacont,1)= mobonface*normcont*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*p(lef)-...
            parameter(1,1,ifacont)*pinterp(facelef1)-parameter(1,2,ifacont)*pinterp(facelef2));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    
end
% Montagem da matriz global

for iface=1:size(inedge,1)
    
    %Define "mobonface" (for "inedge")
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
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    %[mulef, murel]=weightnlfvDMP(kmap,iface);
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxfacelef1=parameter(1,3,ifactual);
    auxfacelef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxfacerel1=parameter(2,3,ifactual);
    auxfacerel2=parameter(2,4,ifactual);
    % calculo dos fluxo parcial a esquerda
    fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*p(lef)-...
                   parameter(1,1,ifactual)*pinterp(auxfacelef1)-parameter(1,2,ifactual)*pinterp(auxfacelef2));
    % calculo dos fluxo parcial a direita        
    fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*p(rel)-...
                       parameter(2,1,ifactual)*pinterp(auxfacerel1)-parameter(2,2,ifactual)*pinterp(auxfacerel2));
    % calculo do fluxo unico na face
    flowrate(iface+size(bedge,1),1)=mobonface*(murel*fluxesq-mulef*fluxdireit);
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
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