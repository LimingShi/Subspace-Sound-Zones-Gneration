function [U_gevd,gevd_values]=Subspace_rotate_2_EVD(RX,RN,U)
RX_NEW=U'*RX*U;RX_NEW=(RX_NEW+RX_NEW')/2;
RN_NEW=U'*RN*U;RN_NEW=(RN_NEW+RN_NEW')/2;
[U_x,gamma_x]=eig_ls_descend(RX_NEW);
inx=find(gamma_x>eps,1,'last');
inx=min(inx,size(U_x,1));
U_x=U_x(:,1:inx(end));
gamma_x=gamma_x(1:inx(end));
sqrt_inv_gamma_x=(1./sqrt(gamma_x));
U_x_scale=U_x.*sqrt_inv_gamma_x';
M=U_x_scale'*RN_NEW*U_x_scale; M=(M+M')/2;
[D,EIG_VALUE_sub_vector]=eig_ls_ascend(M);
inx2=find(EIG_VALUE_sub_vector>eps,1);
z_min=D(:,inx2:end).*sqrt_inv_gamma_x;
U_gevd=U*U_x*z_min;
gevd_values=1./EIG_VALUE_sub_vector(inx2:end);
U_gevd=U_gevd.*(sqrt(gevd_values)');
U_gevd=[U_gevd,zeros(size(U,1),size(U,2)-size(U_gevd,2))];
gevd_values=[gevd_values;eps*ones(size(U,2)-length(gevd_values),1)];
end