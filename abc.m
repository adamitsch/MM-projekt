function abc()
pkg load control;
  
podatki=[10,20,1.5,9.8];
m=podatki(1);
M=podatki(2);
l=podatki(3);
g=podatki(4);

g=g*-1;
   
korak = 0.01;
cas = 10;
tspan = 0:korak:cas;
korakov = cas/korak;

%ZAČETNO STANJE
Y0 = [0; 0; 1.4; 0];
Y = zeros(length(Y0), length(tspan));
Y(:,1)=Y0;
K=[0 0 0 0];

funkcijaa = @(t,x,u) [ x(2); (u+m*l*(sin(x(3))*x(4)^2 ) - m*g*cos(x(3))*sin(x(3)) )/(M+m-m*cos(x(3))^2) ;x(4); (u*cos(x(3))-(M+m)*g*sin(x(3)) + m*l*(cos(x(3))*sin(x(3)) )*x(4)*x(4) )/(m*l*cos(x(3))^2 - (M+m)*l ) ];

for k=2:length(tspan)
   
   %u = -K * Y(:,k-1);
   u = 0;
  
   k1 = korak * funkcijaa(tspan(k), Y(:,k-1), u);
   k2 = korak * funkcijaa(tspan(k) + korak/2, Y(:,k-1) + k1/2, u);
   k3 = korak * funkcijaa(tspan(k) + korak/2, Y(:,k-1) + k2/2, u);
   k4 = korak * funkcijaa(tspan(k) + korak, Y(:,k-1) + k3, u);
   Y(:,k) = Y(:,k-1) + (1/6)*( k1 + 2*k2 + 2*k3 + k4);
   
endfor

%%{
  plot(tspan, Y(1,:) ,'r;pozicijax;' )
  hold on
  plot(tspan, Y(2,:), 'g;hitrost;')
  hold on
  plot(tspan, Y(3,:), 'b;odklon utezi;')
  hold on
  plot(tspan, Y(4,:), 'k;kotna hitrost;')
%}
  
%%{  
for i=1:10:length(Y)
  %
  vektor = Y(:,i);
  x = vektor(1);
  th = vektor(3);
  th = th;
  
  %ali
  %th =vektor(3); 
  
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

  l=1.8;
  
  px = x + l*sin(th);
  py = y - l*cos(th);

  plot([-10 10],[0 0],'g','LineWidth',2)
  hold on
  rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1])
  rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
  rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

  plot([x px],[y py],'w','LineWidth',2)
  rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1])

  % set(gca,'YTick',[])
  % set(gca,'XTick',[])
  xlim([-6 6]);
  ylim([-2 3]);
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