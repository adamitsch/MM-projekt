function aaa()

pkg load control;
  
m = 1;
M = 5;
L = 2;
g = -10;
d = 0;

s = 1;

A = [0 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0];

B = [0; 1/M; 0; s*1/(M*L)];

korak = 0.01;
cas = 20;
tspan = 0:korak:cas;
korakov = cas/korak;

tspan = 0:.001:10;


p = [-1.5; -1.5; -1.6; -1.7];
p = [-1.1; -1.2; -1.3; -1.4];
K = place(A,B,p)


y0 = [0; 0; pi+0.9; 0];

[t,y] = ode45(@(t,y)vozicek(y,m,M,L,g,d,-K*(y-[1; 0; pi; 0])),tspan,y0);

[t y];

%{
  plot(t, y(:,1) ,'r;pozicijax;' )
  hold on
  plot(t, y(:,2), 'g;hitrost;')
  hold on
  plot(t, y(:,3), 'b;odklon utezi;')
  hold on
  plot(t, y(:,4), 'k;kotna hitrost;')
  
  %pause(100);
%}
  
%%{  
Y=y;
for i=1:50:length(t)
  %
  vektor = Y(i,:);
  x = vektor(1);
  th = vektor(3);
  th = th;
 
  %th=mod(-th+pi, 2*pi);
  
  W = 1*sqrt(M/5);  % cart width
  H = .5*sqrt(M/5); % cart height
  wr = .5; % wheel radius
  mr = .5*sqrt(m); % mass radius

  W=2;
  H=1;
  wr=0.5;
  mr=0.7;
  
  % positions
  % y = wr/2; % cart vertical position
  y = wr/2+H/2; % cart vertical position
  w1x = x-.9*W/2;
  w1y = 0;
  w2x = x+.9*W/2-wr;
  w2y = 0;

  l=2;
  
  px = x + l*sin(th);
  py = y - l*cos(th);

  plot([-12 12],[0 0],'g','LineWidth',2)
  hold on
  rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1])
  rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
  rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

  plot([x px],[y py],'w','LineWidth',2)

  rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1])

  % set(gca,'YTick',[])
  % set(gca,'XTick',[])
  xlim([-12 12]);
  ylim([-4 5]);
  set(gca,'Color','k','XColor','w','YColor','w')
  set(gcf,'Position',[10 900 800 400])
  set(gcf,'Color','k')
  set(gcf,'InvertHardcopy','off')   

  % box off
  drawnow
  hold off
  
  if i==1
    pause(1)
  endif

endfor
%}

  
endfunction
  
  
function dy = vozicek(y,m,M,L,g,d,u)

Sy = sin(y(3));
Cy = cos(y(3));
D = m*L*L*(M+m*(1-Cy^2));

dy(1,1) = y(2);
dy(2,1) = (1/D)*(-m^2*L^2*g*Cy*Sy + m*L^2*(m*L*y(4)^2*Sy - d*y(2))) + m*L*L*(1/D)*u;
dy(3,1) = y(4);
dy(4,1) = (1/D)*((m+M)*m*g*L*Sy - m*L*Cy*(m*L*y(4)^2*Sy - d*y(2))) - m*L*Cy*(1/D)*u;

endfunction