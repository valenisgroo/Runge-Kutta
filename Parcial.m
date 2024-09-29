%Valentino Isgró

%Carga de datos
g = load('p1_3k9_01.txt', "-ascii"); 

%Valores iniciales del parcial
b = [0;1];
A = [0,1; -14400, -2.4];
fs=5000;
Dt= 1/fs;
w= 1;
t0=0;
Dim = length(g);
tf = (Dim-1)*Dt;

%Condiciones iniciales del parcial
t = t0:Dt:tf;
tg = zeros(2, Dim);
xg = zeros(2,Dim);
x = zeros(2, Dim);
x(1,1)=0; %x1
x(2,1)=0; %x1

%Runge-Kutta
for j = 1:Dim-1
    k1 = Dt * (A * x(:,j) + b * g(j));
    tg(j) = t(j) + (Dt/(2*w));
    xg(:,j)= x(:,j) + (k1/(2*w));
    k2 = Dt * (A * xg(:,j) + b * g(j));
    x(:,j+1) = x(:,j) + (1-w)*k1 + w*k2;
    if (t(j)==1)
      fprintf('N°de componente j = %.4f\n', j)
      fprintf('x1(2.5) = %.4f\n', x(1,j))
      fprintf('x2(2.5) = %.4f\n', x(2,j))
    endif
endfor

%Gráfica de x1
figure(1)
plot(t, x(1,:),  '-r', 'DisplayName', 'Runge-Kutta')
title(['Gráfico x1 = ', num2str(Dt)])
legend show
grid on


%Derivada por método de la derivada central
x1 = x(1,:);
dx1 = zeros(1, numel(t));
dx1(1)= (1/(2*Dt)) * (-3*x1(1) + 4*x1(2) - x1(3));
dx1(numel(t)) = (1/(2*Dt)) * (3*x1(numel(t)) - 4*x1(numel(t)-1) + x1(numel(t)-2));

for j=2:numel(t)-1
  dx1(j) = (1/(2*Dt)) * (-x1(j-1) + x1(j+1));
endfor

% Integral de trapecios
f11 = x1 .* x1;
IT11 = 0;
for k = 1:Dim-1
  IT11 += (f11(k) + f11(k+1)) * Dt / 2;
endfor
fprintf('Resultado de IT11: %.4f\n', IT11);

% Integral de trapecios
f22 = dx1.*dx1;
IT22 = 0;
for k = 1:Dim-1
  IT22 += (f22(k) + f22(k+1)) * Dt / 2;
endfor
fprintf('Resultado IT22: %.4f\n', IT22);