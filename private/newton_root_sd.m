
function x0=newton_root_sd(f,fprime,x0,maxIterations,tolerance)
if nargin<=3
    maxIterations=1e4;
    tolerance=1e-7;
end

for i = 1 : maxIterations
    
    y = f(x0);
%     if (y>=0)
%         break;
%     end
    
    
    yprime = fprime(x0);
    
    if(abs(yprime) < eps) %Don't want to divide by too small of a number
        % denominator is too small
        break; %Leave the loop
    end
    
    x1=-1; mu=1;
    while x1<0        
        x1 = x0 - mu*y/yprime; %Do Newton's computation
        mu=mu/2;
    end
    
    
    if(norm(x1 - x0)/norm(x0) <= tolerance) %If the result is within the desired tolerance
        break; %Done, so leave the loop
    end
    
    x0 = x1; %Update x0 to start the process again
    
end
end