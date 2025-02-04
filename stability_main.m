disp('####################################################################### ');
disp('##                                                                   ## ');
disp('##                          Chandan Kr                               ## ');
disp('##                   The Ohio State University    	                 ## ');
disp('##     -------------------------------------------------------       ## ');
disp('##                                                                   ## ');
disp('##            Orr Sommerfeld stability - Temporal stability          ## ');
disp('##                Tollmien Schlichting wave - Blasius BL             ## ');
disp('##     -------------------------------------------------------       ## ');
disp('####################################################################### ');

clc; clear; close all;
N = 201;
[u,ddu,y,dy] = blasiusBL(N);

figure(1)
plot(u,y,'-k',LineWidth=2)
xlabel('u')
ylabel('y')
title('Mean boundary layer')
set(gca,'FontSize',14, 'FontWeight','bold')
set(gcf,'Position',[100 100 400 500])

alp =0.179;
R =580;

[c,v]=OS_FDM(u,ddu,dy,R,alp,N);

% plot eigevalue
cr = real(c); ci = imag(c);

figure(2)
  plot(cr,ci,'ok', LineWidth=2)
  hold on
  yline(0,'--r', LineWidth=1)
  xlabel('c_r'); ylabel('c_i');
  ylim([-1 1]); xlim([0 1])
  tt = strcat(sprintf('Eigenvalue spectrum at alpha =%1.3f; Re = %1.1f',alp,R) );
  title(tt,'Interpreter','tex');
  set(gca,'FontSize',14, 'FontWeight','bold')

vp = (flipud([0;0; v(:,1); 0;0]))' ;

figure(3)
  plot(abs(vp),y,'-k','LineWidth',2.0, DisplayName='abs');
  hold on
  plot(imag(vp),y,'-r','LineWidth',2.0,  DisplayName='imag');
  plot(real(vp),y,'-b', 'LineWidth',2.0,  DisplayName='real');
  hold off
  xlabel('$\hat{v}$','Interpreter','latex')
  ylabel('y')
  legend
  title('v Perturbation')
  set(gca,'FontSize',14, 'FontWeight','bold')
  set(gcf,'Position',[100 100 400 500])

  for i=1:N
      if(i==1)
          dvp(i) = (vp(i+1)-vp(i))/dy;  
      elseif(i==N)
          dvp(i) = (vp(i)-vp(i-1))/dy; 
      else
          dvp(i) = (1/2*vp(i+1)-1/2*vp(i-1))/dy;
      end
  end
 
  up =-dvp./(1i*alp);

figure(4)
  plot(abs(up),y,'-k','LineWidth',2.0, DisplayName='abs');
  hold on
  plot(imag(up),y,'-r','LineWidth',2.0,  DisplayName='imag');
  plot(real(up),y,'-b', 'LineWidth',2.0,  DisplayName='real');
  hold off
  xlabel('$\hat{u}$','Interpreter','latex')
  ylabel('y')
  legend
  title('u Perturbation')
  set(gca,'FontSize',14, 'FontWeight','bold')
  set(gcf,'Position',[100 100 400 500])

% filter eigenvalue and plot corresponding eigenvector:
% Neutral curve 

alp =0.01:0.0025:0.21;
R = 300:100:8000;

for i = 1:length(alp)
    for j =1:length(R)
        [c,v]=OS_FDM(u,ddu,dy,R(j),alp(i),N);

        dc=find(real(c) > 0 & real(c) < 1 & imag(c) > -1 & imag(c) < 1);
        if ~isempty(dc)
          index=find(imag(c)==max(imag(c(dc))));
          ciMin(i,j)=c(index);
        else
          ciMin(i,j)=NaN;
        end
        clear c

        clc
        iter =  (((i*length(R)) - length(R)+j)/(length(alp)*length(R)))*100;
        fprintf('\n alpha=%f. Progress %1.2f%%... (Maximum alpha =%f)\n', alp(i), iter, alp(end))

     end
end
   
figure(5) 
  contourf(R,alp,imag(ciMin));
  xlabel('Re'); ylabel('\alpha')
  ylim([0 0.25]); xlim([0 8000])
  title('Neutral curve')
  set(gca,'FontSize',14, 'FontWeight','bold')
  colorbar;


