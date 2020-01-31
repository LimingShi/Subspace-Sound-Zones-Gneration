function [x,R_total,diag_values,alpha_set,beta_set] = conjgrad(A, b, x,iteration_num)
N=length(b);
R_total=zeros(N,iteration_num);
r = b - A * x;
p = r;
rsold = r' * r;
R_total(:,1)=p;


diag_values=zeros(iteration_num,1);
alpha_set=zeros(iteration_num,1);
beta_set=zeros(iteration_num,1);

for i = 1:iteration_num
    Ap = A * p;
    diag_values(i,1)=p'*Ap;
    
    
    alpha = rsold / diag_values(i,1);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    %     if sqrt(rsnew) < 1e-10
    %         break;
    %     end
    beta_set(i)=rsnew / rsold ;
    
    p = r + beta_set(i) * p;
    rsold = rsnew;
    R_total(:,i+1)=p;
    
    
    alpha_set(i,1)=alpha;
end

end