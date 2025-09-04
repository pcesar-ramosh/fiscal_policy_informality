function [SupI, DemI, SupF, DemF, SupTot, DemTot, Stot] = ...
    market_curves_by_r_BASE(share_I, sI_vec, sF_vec, trans, cfg, rgrid)

K = numel(rgrid);
SupI=zeros(K,1); DemI=zeros(K,1);
SupF=zeros(K,1); DemF=zeros(K,1);
SupTot=zeros(K,1); DemTot=zeros(K,1); Stot=zeros(K,1);

for k=1:K
    r_ = rgrid(k);
    [~, ~, g_all, ~, a_g] = huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg, r_);
    da = a_g(2)-a_g(1);
    supI=0; demI=0; supF=0; demF=0;
    for j=1:numel(sI_vec)
        g = g_all{j};
        supI = supI + sum(max(a_g,0).*g(:,1))*da;
        demI = demI + sum(-min(a_g,0).*g(:,1))*da;
        supF = supF + sum(max(a_g,0).*g(:,2))*da;
        demF = demF + sum(-min(a_g,0).*g(:,2))*da;
    end
    SupI(k)=supI; DemI(k)=demI; SupF(k)=supF; DemF(k)=demF;
    SupTot(k)=supI+supF; DemTot(k)=demI+demF; Stot(k)=SupTot(k)-DemTot(k);
end
end
