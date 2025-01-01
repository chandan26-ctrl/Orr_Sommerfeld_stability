



%  D   B    A    C   E   N-1    I =1
%  B    A    C   E       N-2   I= 2
 % D   B    A    C   E         I = 3
           ..................
                       %  D   B    A    C   E        I = x-2
                           %  D   B    A    C   3     I = x-1
                        % D   B    A    C   E   2     I =x

% J+2  J+1  J    J-1  J-2



function [c,v]=OS_FDM(u,ddu,dy,R,alp,N)
RE =R;
for j=3:N-2
    E(j) = 1i/(RE*alp)/dy^4;   
    C(j) = 1i/(RE*alp)*(-4/dy^4 - 2*alp^2/dy^2) +u(j)/dy^2;
    A(j) = 1i/(RE*alp)*(6/dy^4 + 4*alp^2/dy^2 + alp^4) -2*u(j)/dy^2- u(j)*alp^2-ddu(j);
    B(j) = 1i/(RE*alp)*(-4/dy^4 - 2*alp^2/dy^2) +u(j)/dy^2;
    D(j) = 1i/(RE*alp)/dy^4;
end

E = E(3:end); C = C(3:end); A = A(3:end); B = B(3:end); D = D(3:end);

% j = N-1;
%     E(j) = 1i/(RE*alp)/dy^4;   
%     C(j) = 1i/(RE*alp)*(-4/dy^4);
%     A(j)  = 1i/(RE*alp)*(6/dy^4);
%     B(j) = 1i/(RE*alp)*(-4/dy^4 - 2*alp^2/dy^2) +u(j)/dy^2;
%     D(j) = 1i/(RE*alp)*(1/dy^4 + 4*alp^2/dy^2+alp^4) + u(j)*(-2/dy^2-alp^2) -ddu(j);
% 
% j =2;
%     E(j) = 1i/(RE*alp)*(1/dy^4 + 4*alp^2/dy^2+alp^4) + u(j)*(-2/dy^2-alp^2) -ddu(j);   
%     C(j) = 1i/(RE*alp)*(-4/dy^4 - 2*alp^2/dy^2) +u(j)/dy^2;
%     A(j) = 1i/(RE*alp)*(6/dy^4);
%     B(j) = 1i/(RE*alp)*(-4/dy^4);
%     D(j) = 1i/(RE*alp)/dy^4; 


   % N =100

    sz = size(E); N1 = sz(1,2);

    AMT = zeros(N1,N1);

    for i =1:N1
         AMT(i,i) = A(N1+1-i);
    end
     k=1;
    for i =1:N1-1
         k= k+1;
        AMT(i,k)= C(N1+1-i);
    end

    k=2;
    for i =1:N1-2
         k= k+1;
        AMT(i,k)= E(N1+1-i);
    end

    k=0;
    for i =2:N1
         k= k+1;
        AMT(i,k)= B(N1+1-i);
    end

    k=0;
    for i =3:N1
         k= k+1;
        AMT(i,k)= D(N1+1-i);
    end
    
    % M O      N-1  1
    % L M O    N-2  2


         % L M 0   3  X-1
              % L M  2   X
  for j=3:N-2
      BCAP(j) = 1/dy^2;
      ACAP(j) = -2/dy^2-alp^2;
      CCAP(j) = 1/dy^2;
  end

  BCAP = BCAP(3:end); ACAP = ACAP(3:end); CCAP = CCAP(3:end); 


   BMT = zeros(N1,N1);

    for i =1:N1
         BMT(i,i) = ACAP(N1+1-i);
    end

    k=1;
    for i =1:N1-1
         k= k+1;
        BMT(i,k)= CCAP(N1+1-i);
    end

        k=0;
    for i =2:N1
         k= k+1;
        BMT(i,k)= BCAP(N1+1-i);
    end



[v,ei]=eig(AMT, BMT); % solves eigenvalue problem
c =diag(ei);
end
