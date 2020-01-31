function [eveMtx_eig,evaVec_eig]=eig_ls_ascend(X)
[eveMtx_eig, evaVec_eig] = eig(X);
evaVec_eig=diag(evaVec_eig);

[evaVec_eig,inx]=sort(evaVec_eig,'ascend');

eveMtx_eig=eveMtx_eig(:,inx);

end