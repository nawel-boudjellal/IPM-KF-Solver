%The step size%
function pasp=depl(x,deltax)
%min=0.99;
min=1;
n=length(x);
for i=1:n
    if deltax(i)<0
        if min>(-x(i)/deltax(i))
            min=-x(i)/deltax(i);
        end
    end
end
pasp=min;
end
        
