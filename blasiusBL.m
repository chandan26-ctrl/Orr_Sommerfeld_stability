function [u,ddu,y,h] = blasiusBL(N)
L = 20;
h = L/(N-1);
%% Numerical Solution of Blasius Equation Using Runge-Kutta
f1 = @(x, y1, y2, y3) y2;
f2 = @(x, y1, y2, y3) y3;
f3 = @(x, y1, y2, y3) -0.5*y1*y3;
eta = 0:h:L;
y1(1) = 0;
y2(1) = 0;
y3(1) = 0.33206;    %f''(0) =  0.33206 (Boundary Layer Theory, H. Schlichting) 
for i = 1:(length(eta)-1)
a = h.*[f1(eta(i), y1(i), y2(i), y3(i)), f2(eta(i), y1(i), y2(i), y3(i)), f3(eta(i), y1(i), y2(i), y3(i))];
b = h.*[f1(eta(i), y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f2(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f3(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2)];
c = h.*[f1(eta(i), y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f2(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f3(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2)];
d = h.*[f1(eta(i), y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f2(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f3(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3))];
y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end
ddu = - (y1.*y3)/2;
u = y2;
y = eta;
end

