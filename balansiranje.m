function balansiranje(abc)
 pkg load control;
 
 
%{ 
 KOMENTAR
 m - masa uteži
 M - masa vozička
 l - dolžina palice
 g - težni pospešek
%}

podatki=[10,5,1,9.8];
%podatki=[20,100,1,9.8];
%podatki=[50,10,1,9.8];

m=podatki(1);
M=podatki(2);
l=podatki(3);
g=podatki(4);

% J, b -- J je vztrajnostni moment, b trenje
% J =  ** kg m^2 / s^2
% b =  ** Ns / m
b=0.1;

J=100;
J = (1/3) * m * l * l

% J = 1/3 * m * l^2


%!!
%{
M=500;
m=1;
J=100;
l=10;
g=9.8;
K=[-1 120 -4 200];
%}

p = J * (M+m) + M*m*l^2; % denominator ?

% -------  matriki A in B --------------

% x theta x' thet'a ----> x x' theta theta'

%%{
A = [ 0         1               0       0;
      0 (-(J+m*l^2)*b)/p  (m^2*g*l^2)/p  0;
      0         0               0       1;
      0    (-m*l*b)/p   (m*g*l*(M+m))/p 0];

A_prava=A;
A_leva = A;
A_leva(4,2)=-1*A_leva(4,2);
A_leva(4,3)=-1*A_leva(4,3);



B = [   0;
   (J+m*l^2)/p;
        0;
   (m*l)/p];
   
  lastneA= eig(A)
  
%}

% -------------

%%{   
A = [ 0 1           0           0;
      0 0       (-g*m)/M        0;
      0 0           0           1;
      0 0 (g*(l+M/m))/(l*(M/m)) 0];

% *neka fora spodaj/zgoraj!      
%Adruga = A;
%Adruga = -Adruga(4,3);
%A(4,3)=0;

B = [   0;
      1/M;
        0;
   1/(M*l)];

C = [1 0 0 0];   
%}

korak = 0.01;
cas = 10;
tspan = 0:korak:cas;  %čas
korakov = cas/korak;

% x = [ pozicija, hitrost odklon uteži kotna hitrost]'
Y0 = [2;0.5;0.1;0.01];
Y0 = [5;0;2;0];
Y0 = [1;0;0.5;0];

% x theta x' theta'  za testiranje po funkciji iz pdfja 
Y0 = [1; 0.5 ; 0 ; 0];

%Y0 = [0 1 0 0];
Y0 = [0; 0; 0.6; 0];
Y0 = [0; 0; 2; 0];
%Y0 = [1 0 1 0];

%Y0 = [1,0,pi/2,0];
%Y0=abc

Y = zeros(length(Y0), length(tspan));
Y(:,1)=Y0;

% ======= vektor K  ======================
% u(t) = -K * x(t) = - (K1x1 + K2x2 + K3x3 + K4x4)
%K = [ -10 -20 -30 -40];
%K = [ -1 -5 100 200];

% iz lastnih vrednosti - rečeš kakšne lastne vrednosti hočeš -> funkcija place
p = [-1.2; -1.2; -1.2; -1.2];
%p = [-0.5; -0.5; -0.5; -0.5];
%mal bl agresivno
%p = [-2,-2,-2,-2];

Q=[1 0 0 0; 0 1 0 0; 0 0 10 0; 0 0 0 100];
R = 0.01;

%kontroliranje

%K = [ -1 -5 100 200];
%K = [ -1 -4 120 20];
%K = [-15 -50 -2000 500];
K = place(A,B,p);
%K = lqr(A,B,Q,R);
%K = [0 0 0 0];


% =========== RUNGE KUTTA ============================


%funkcija = @(t,x)  ([x(2);x(4);inv([M+m m*l*cos(x(3)); J+m*l^2 m*l*cos(x(3))])*[-b*x(2)+m*l*sin(x(3))*x(4)^2+K*x; -m*g*l*sin( x(3) ) ] ]' * [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])';

funkcija = @(t,x,K) [x(2);  x(4);  inv([M+m m*l*cos(x(3)); J+m*l^2 m*l*cos(x(3))])  *  [-b*x(2)+m*l*sin(x(3))*x(4)^2+K*x;   -m*g*l*sin( x(3) ) ] ];

funkcijaa = @(t,x,u) [x(2); (u+m*l*(sin(x(3))*x(4)^2 ) - m*g*cos(x(3))*sin(x(3)) )/(M+m-m*cos(x(3))^2) ;x(4); (u*cos(x(3))-(M+m)*g*sin(x(3)) + m*l*(cos(x(3))*sin(x(3)) )*x(4)*x(4) )/(m*l*cos(x(3))^2 - (M+m)*l ) ];



f = @(t,x,u,A) A*x + B*u;


% pr funkcijah obrnemo tud K



%K = [ -4.0816   1127.4422    -16.4265    181.7687];
K = [0,0,0,0];



for k=2:length(tspan)

  %katero matriko uporabit
  %kot=Y(3,k-1);
  %kot=mod(kot,2*pi);
  
  %u = 0; 
  u = -K * Y(:,k-1);

%================ matriki A & B - RK3 ============ 
  %{
  k1 = A*Y(:,k-1) + B*u; 
  q1 = Y(:,k-1) + k1*korak/2;
  
  k2 = A*q1 + B*u;
  Y(:,k) = Y(:,k-1) + k2*korak;
  
  %rk http://lpsa.swarthmore.edu/NumInt/NumIntSecond.html  
  
  %}
%================= iz pdfja ================
   
   
   %{
   k1 = korak * f(tspan(k), Y(:,k-1), u ,A);
   k2 = korak * f(tspan(k) + korak/2, Y(:,k-1) + k1/2 , u,A);
   k3 = korak * f(tspan(k) + korak/2, Y(:,k-1) + k2/2, u,A);
   k4 = korak * f(tspan(k) + korak, Y(:,k-1) + k3 ,u ,A);
   
   Y(:,k) = Y(:,k-1) + (1/6)*( k1 + 2*k2 + 2*k3 + k4);
   %}
   
   %{
   K = K * [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
   u = -K * Y(:,k-1);
   
   k1 = korak * funkcija(tspan(k), Y(:,k-1), K);
   k2 = korak * funkcija(tspan(k) + korak/2, Y(:,k-1) + k1/2, K);
   k3 = korak * funkcija(tspan(k) + korak/2, Y(:,k-1) + k2/2, K);
   k4 = korak * funkcija(tspan(k) + korak, Y(:,k-1) + k3, K);
   Y(:,k) = Y(:,k-1) + (1/6)*( k1 + 2*k2 + 2*k3 + k4);
   
   %}
   
   %{
   k1 = korak * funkcijaa(tspan(k), Y(:,k-1), u);
   k2 = korak * funkcijaa(tspan(k) + korak/2, Y(:,k-1) + k1/2, u);
   k3 = korak * funkcijaa(tspan(k) + korak/2, Y(:,k-1) + k2/2, u);
   k4 = korak * funkcijaa(tspan(k) + korak, Y(:,k-1) + k3, u);
   Y(:,k) = Y(:,k-1) + (1/6)*( k1 + 2*k2 + 2*k3 + k4);
   %}
  
   
  
endfor

%========== runge kutta plot ===========
% za @f tazgornji, za @funkcija taspodnji plot

%{
plot(tspan, Y(1,:) ,'r;pozicijax;' )
hold on
plot(tspan, Y(2,:), 'g;hitrost;')
hold on
plot(tspan, Y(3,:), 'b;odklon utezi;')
hold on
plot(tspan, Y(4,:), 'k;kotna hitrost;')
%}

% ta je za funkcijo
%{
plot(tspan, Y(1,:) ,'r;pozicijax;' )
hold on
plot(tspan, Y(2,:), 'b;odklon;')
hold on
plot(tspan, Y(3,:), 'g;hitrost;')
hold on
plot(tspan, Y(4,:), 'k;kotna hitrost;')
%}

% ============= ode45  & plot==========

%%{
K = [0,0,0,0];
K = place(A,B,p);

funkcija = @(t,x)  ([x(2);x(4);inv([M+m m*l*cos(x(3)); J+m*l^2 m*l*cos(x(3))])*[-b*x(2)+m*l*sin(x(3))*x(4)^2+K*x; -m*g*l*sin( x(3) ) ] ]' * [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])';
t0 = 0:.01:10;
[T,Y] = ode45(funkcija, [tspan], [0; 0; pi; 0]);
numel(T);
numel(Y);
%y0 = interp1(T,Y,t0);

%funkcijaMatrike = @(t,x) ( A * x + B * (K*x));
%[T,Y] = ode45(funkcijaMatrike, [0,10], [1; 0; 0.5; 0]);

plot(T, Y(:,1) ,'r;pozicijax;' );
hold on
plot(T, Y(:,2), 'b;odklon utezi;');
hold on
plot(T, Y(:,3), 'g;hitrost;');
hold on
plot(T, Y(:,4), 'k;kotna hitrost;');
%}
 


% ============= RUNGE KUTTA iz učilnce ========
%{
function [t, Y] = rk4(f, interval, Y0, h)
%[t, Y] = rk4(f, [t0, tk], Y0, h) resi DE
%Y' = f(Y, t) pri zacetnem pogoju Y(t0) = Y0 s standardno 
%Runge-Kutta metodo 4. reda s korakom h na intervalu [t0, tk].

%pripravimo vrstico casov t
t = interval(1):h:interval(2);

%pripravimo vrednosti Y
Y = zeros(length(Y0), length(t));
Y(:, 1) = Y0;

%pozenemo Runge-Kutta metodo
for k = 2:length(t)
	%poracunamo vrednosti k1, ..., k4 in...
	k1 = h*f(Y(:, k-1), t(k-1));
	k2 = h*f( Y(:, k-1) + k1/2, t(k-1) + h/2);
	k3 = h*f(Y(:, k-1) + k2/2, t(k-1) + h/2);
	k4 = h*f(Y(:, k-1) + k3, t(k-1) + h);
	%... dodamo utezeno povprecje teh vrednosti Y
	Y(:, k) = Y(:, k-1) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
end

%}


Y(:,1)
Y(:,2)
%%{
for i=1:5:length(Y)
  %
  vektor = Y(:,i);
  vektor = Y(i,:);
  
  
  x = vektor(1);
  th = vektor(3);
  th = vektor(2);
  th = th+pi;
  
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

  l=1.5;
  
  px = x + l*sin(th);
  py = y - l*cos(th);

  plot([-10 10],[0 0],'w','LineWidth',2)
  hold on
  rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1])
  rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
  rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

  plot([x px],[y py],'w','LineWidth',2)

  rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1])

  % set(gca,'YTick',[])
  % set(gca,'XTick',[])
  xlim([-10 10]);
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
