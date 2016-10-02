function [flowrate,flowresult]=FflowrateNLFV(p, pinterp, parameter,mobility)
global inedge coord bedge bcflag centelem phasekey numcase
if numcase==31.2
        auxbcflag=bcflag;
        auxbcflag(2,2)=-1;
    else
        auxbcflag=bcflag;
 end
bedgesize=size(bedge,1);
inedgesize=size(inedge,1);
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
        flowrate(ifacont,1)= normcont*auxbcflag(r,2);
    else

        flowrate(ifacont,1)=mobonface*normcont*(parameter(1,1,ifacont)*(p(lef)-pinterp(parameter(1,3,ifacont)))+...
                                                    parameter(1,2,ifacont)*(p(lef)-pinterp(parameter(1,4,ifacont))));        
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);  
    
end

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
    
    %% calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    alef=norma*(parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual)));
    % direita
    
    arel= norma*(parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual)));
    %% calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    if alef==0 && arel==0
        mulef= 0.5;
        murel=1-mulef;
    else
        mulef=abs(arel)/(abs(alef)+abs(arel));
        murel=1-mulef;
    end
    %% calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ALR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    
   flowrate(iface+size(bedge,1),1)=mobonface*(ALL*p(lef)-ALR*p(rel)); 
   
   %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);  
    
end

end