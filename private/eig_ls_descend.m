function [eveMtx_eig,evaVec_eig]=eig_ls_descend(X)
[eveMtx_eig, evaVec_eig] = eig(X);
evaVec_eig=diag(evaVec_eig);

[evaVec_eig,inx]=sort(evaVec_eig,'descend');

eveMtx_eig=eveMtx_eig(:,inx);

end