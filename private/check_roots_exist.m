function mu_final=check_roots_exist(c_square_1,c_square_2,EIG_VALUE_total,K1,K2)

%% for the second source; starting point
mu_vaule11_try_2=0;
a_value=sum(c_square_2.*EIG_VALUE_total{2}./((EIG_VALUE_total{2}+mu_vaule11_try_2).^2));
f_try=@(mu_vaule2) a_value/K2-sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule2)./(EIG_VALUE_total{1}+mu_vaule2));
if f_try(0)<0
    f_try_prime=@(mu_vaule2) 2*sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule2)./(EIG_VALUE_total{1}+mu_vaule2)./(EIG_VALUE_total{1}+mu_vaule2));
    x0=0;
    initvalue_u11=newton_root_sd(f_try,f_try_prime,x0);
    initvalue_u2=0;
    
else
    initvalue_u11=0;
    b_value=sum(c_square_1./(EIG_VALUE_total{1}+initvalue_u11)./(EIG_VALUE_total{1}+initvalue_u11));
    f_try2=@(mu_vaule2) b_value*K2-sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    f_try2_prime=@(mu_vaule2) 2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    x0=0;
    initvalue_u2=newton_root_sd(f_try2,f_try2_prime,x0);
    
%     u1u2=[initvalue_u11,initvalue_u2]
end
u1u2_source2=[initvalue_u11,initvalue_u2];



%% for the first source starting point
mu_vaule11_try_2=u1u2_source2(1);
a_value=sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule11_try_2)./(EIG_VALUE_total{1}+mu_vaule11_try_2));
f_try=@(mu_vaule2) a_value/K1-sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));

if f_try(0)<0
    f_try_prime=@(mu_vaule2) 2*sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    x0=0;
    initvalue_u1=newton_root_sd(f_try,f_try_prime,x0);
else
    initvalue_u1=0;
end
u1u2_source1=[u1u2_source2(1),initvalue_u1];


if initvalue_u2>=initvalue_u1;
    mu_final=[u1u2_source2(1),initvalue_u1];
else
    mu_vaule11_try=1e4;
    a_value=sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule11_try)./(EIG_VALUE_total{1}+mu_vaule11_try));
    f_try=@(mu_vaule2) a_value/K1-sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    f_try_prime=@(mu_vaule2) 2*sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    x0=0;
    finalv__=newton_root_sd(f_try,f_try_prime,x0);
    
    
    
    
    b_value=sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule11_try)./(EIG_VALUE_total{1}+mu_vaule11_try));
    f_try2=@(mu_vaule2) b_value*K2-sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    f_try2_prime=@(mu_vaule2) 2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
    x0=0;
    finalv2____=newton_root_sd(f_try2,f_try2_prime,x0);
    
    if  finalv__<finalv2____ 
        
        f1=@(mu_vaule1,mu_vaule2) sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))-K1;
        f2=@(mu_vaule1,mu_vaule2) sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))-K2;
        
        
        g11=@(mu_vaule1,mu_vaule2) -2*sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2));
        g12=@(mu_vaule1,mu_vaule2) 2*sum(c_square_1.*EIG_VALUE_total{1}./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))...
            *sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/(sum(c_square_2./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)))^2;
        g21=@(mu_vaule1,mu_vaule2) 2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))...
            *sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1))/(sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1)))^2;
        g22=@(mu_vaule1,mu_vaule2) -2*sum(c_square_2.*EIG_VALUE_total{2}./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2)./(EIG_VALUE_total{2}+mu_vaule2))/sum(c_square_1./(EIG_VALUE_total{1}+mu_vaule1)./(EIG_VALUE_total{1}+mu_vaule1));
        x0=[0;0];
        
        mu_final=newton_root_sd_dim2(f1,f2,g11,g12,g21,g22,x0);
    else
        disp('No solution');
        mu_final=[nan,nan];
        
    end
    
end



% 
% if finalv__<finalv2____ && initvalue_u2(1) <=initvalue_u1
%     RESULTS=1;
% else
%      RESULTS=0;
% %     disp('solution exists')
% % else
% %     disp('solution do not exist')
% end
end