function balansiranje(abc)
%{ 
 KOMENTAR
 m - masa uteži
 M - masa vozička
 l - dolžina palice
 g - težni pospešek
%}

podatki=[0.1,1,1,9.81];
%podatki=[5,100,1,9.8];
%podatki=[80,10,1,9.8];

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

%{   
A = [ 0 1           0           0;
      0 0       (-g*m)/M        0;
      0 0           0           1;
      0 0 (g*(l+M/m))/(l*(M/m)) 0]

% *neka fora spodaj/zgoraj!      
%Adruga = A;
%Adruga = -Adruga(4,3);
%A(4,3)=0;

B = [   0;
      1/M;
        0;
   1/(M*l)]

C = [1 0 0 0];   
%}

korak = 0.01;
cas = 5;
tspan = 0:korak:cas;  %čas
korakov = cas/korak;

% x = [ pozicija, hitrost odklon uteži kotna hitrost]'
Y0 = [2;0.5;0.1;0.01];
Y0 = [5;0;2;0];
Y0 = [1;0;0.5;0];

% x theta x' theta'  za testiranje po funkciji iz pdfja 
Y0 = [1; 0 ; 0.5 ; 0];

%Y0 = [1,0,pi/2,0];
%Y0=abc

Y = zeros(length(Y0), length(tspan));
Y(:,1)=Y0;

% ======= vektor K  ======================
% u(t) = -K * x(t) = - (K1x1 + K2x2 + K3x3 + K4x4)
K = [ -10 -20 -30 -40];
K = [ -1 -5 100 200];

% iz lastnih vrednosti - rečeš kakšne lastne vrednosti hočeš -> funkcija place
p = [-1; -1; -1; -1];
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

funkcija = @(t,x) [x(2);  x(4);  inv([M+m m*l*cos(x(3)); J+m*l^2 m*l*cos(x(3))])  *  [-b*x(2)+m*l*sin(x(3))*x(4)^2+K*x;   -m*g*l*sin( x(3) ) ] ];

f = @(t,x,u,A) A*x + B*u;


% pr funkcijah obrnemo tud K

K = K * [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1]

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
   
   %%{
   
  
   k1 = korak * funkcija(tspan(k), Y(:,k-1));
   k2 = korak * funkcija(tspan(k) + korak/2, Y(:,k-1) + k1/2);
   k3 = korak * funkcija(tspan(k) + korak/2, Y(:,k-1) + k2/2);
   k4 = korak * funkcija(tspan(k) + korak, Y(:,k-1) + k3);
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

%{
K = [0,0,0,0];


funkcija = @(t,x)  ([x(2);x(4);inv([M+m m*l*cos(x(3)); J+m*l^2 m*l*cos(x(3))])*[-b*x(2)+m*l*sin(x(3))*x(4)^2+K*x; -m*g*l*sin( x(3) ) ] ]' * [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])';

[T,Y] = ode45(funkcija, [0,5], [1; 0; 0.5; 0]);

%funkcijaMatrike = @(t,x) ( A * x + B * (K*x));
%[T,Y] = ode45(funkcijaMatrike, [0,10], [1; 0; 0.5; 0]);

plot(T, Y(:,1) ,'r;pozicijax;' );
hold on
plot(T, Y(:,2), 'g;hitrost;');
hold on
plot(T, Y(:,3), 'b;odklon utezi;');
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

endfunction
