clc
clear
close all


%% Parámetros del sistema
a = 0.2; b = 0.2; c = 5.7;

% Condiciones iniciales --> busco un punto donde pueda ver el
% comportamiento caotico
xi_0 = [0.01, -0.03, 0.04]; 
xi_1 = [6, -28, 28];
xi_2 = [1 1 1]; 


% Ventana temporal
t_span = linspace(0, 200, 5000);

% Sistema de Rössler

rossler = @(t,x) [-x(2) - x(3);
           x(1) + a*x(2);
           b + x(3)*(x(1) - c)];

% Resolución
[t, x] = ode45(rossler, t_span, xi_0);
[~, y] = ode45(rossler, t_span, xi_1);
[~, z] = ode45(rossler, t_span, xi_2);

%% Diagrama de fase
figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(y(:, 1), y(:, 2), y(:, 3), 'b--', 'LineWidth', 1.5);
plot3(z(:, 1), z(:, 2), z(:, 3), 'r-', 'LineWidth', 1.5);
plot3(xi_0(1), xi_0(2), xi_0(3), 'r*', 'MarkerSize',10,'MarkerFaceColor', 'g'); % Inicio X01
plot3(xi_1(1), xi_1(2), xi_1(3), 'go', 'MarkerFaceColor', 'b'); % Inicio x1
plot3(xi_2(1), xi_2(2), xi_2(3), 'bo', 'MarkerFaceColor', 'r'); % Inicio x2
xlabel('x'); ylabel('y'); zlabel('z');
title('Diagrama de Fase del Sistema No Lineal');
legend('Trayectoria 1', 'Trayectoria 2', 'Trayectoria 3', 'Inicio 1', 'Inicio 2', 'Inicio 3');
grid on; hold off;

% Gráfico de trayectorias temporales sistema no linealizado
figure(2)
subplot(3,1,1)
plot(t, x(:,1), 'g'); 
xlabel('Tiempo');
ylabel('x(t)');
xlim([0 200]);
title('Trayectorias Temporales de las soluciones ODE - No Lineal');
grid on
subplot(3,1,2)
plot(t, x(:,2), 'b'); %
xlabel('Tiempo');
ylabel('y(t)');
xlim([0 200]);
grid on
subplot(3,1,3)
plot(t, x(:,3), 'r'); 
ylabel('z(t)');
xlim([0 200]);
grid on;

%% Puntos de equilibrio
F = @(x) [-x(2) - x(3);
           x(1) + a*x(2);
           b + x(3)*(x(1) - c)];

f=@(z) a*z^2-c*z+b; %3ra ecuacion despejada

z_1=fzero(f,0);z_2=fzero(f,20);

xeq1= [z_1*a;-z_1;z_1]; xeq2= [z_2*a;-z_2;z_2];

disp('Puntos de equilibrio:')
disp(xeq1')
disp(xeq2')

%% Matriz Jacobiana

syms x y z

f1 = -y - z;
f2 = x + a*y;
f3 = b + z*(x - c);

Fv = [f1; f2; f3];

% Jacobiana simbólica
J = jacobian(Fv, [x y z]);
disp('Matriz Jacobiana simbólica:')
disp(J)

%% Linearizacion del sistema
disp('Jacobiana en el punto de equilibrio x1:')
J_eval1 = double(subs(J, [x y z], xeq1'));
disp(J_eval1)

f_lin1 = @(t, X) J_eval1 * X;

disp('Jacobiana en el punto de equilibrio x2:')
J_eval2 = double(subs(J, [x y z], xeq2'));
disp(J_eval2)

f_lin2 = @(t, X) J_eval2 * X;


t_span = linspace(0, 50, 2000); %para sistema lineal necesito tiempo más corto

% Simulación
[t, X1] = ode45(f_lin1, t_span, xeq1);
[~, X2] = ode45(f_lin2, t_span, xeq2);

% Gráfico
figure(3);
subplot(2,1,1)
plot(t, X1);
xlabel('Tiempo'); ylabel('Estados');
title('Comportamiento del sistema linealizado alrededor de x1');
legend('x', 'y', 'z');
grid on;
subplot(2,1,2)
plot(t, X2);
xlabel('Tiempo'); ylabel('Estados');
title('Comportamiento del sistema linealizado alrededor de x2');
legend('x', 'y', 'z');
grid on;

%% Autovalores
lambda1 = eig(J_eval1);
lambda2 = eig(J_eval2);
disp('Autovalores calculados en MATLAB:')
disp('Para el primer punto de equilibrio')
disp(lambda1)
disp('Para el segundo punto de equilibrio')
disp(lambda2)

%% Simulaciones para modelo linealizado - Diagramas de Fase

% Repetición simulaciones con condiciones iniciales originales

% Condiciones iniciales originales (columna para restar)
xi_0_col = xi_0';  
xi_1_col = xi_1';
xi_2_col = xi_2';

% Llevar condiciones iniciales al sistema de coordenadas de xeq1 - mover
% centro de referencia

% (perturbación = C.I. originales - punto de equilibrio)
X0_lin1_a = xi_0_col - xeq1; % C.I. 1 (desde xi_0)
X0_lin1_b = xi_1_col - xeq1; % C.I. 2 (desde xi_1)
X0_lin1_c = xi_2_col - xeq1; % C.I. 3 (desde xi_2)

% Simulación de las 3 trayectorias para diferentes CI, para el sistema linealizado f_lin1
[t1a, X1a] = ode45(f_lin1, t_span, X0_lin1_a);
[t1b, X1b] = ode45(f_lin1, t_span, X0_lin1_b);
[t1c, X1c] = ode45(f_lin1, t_span, X0_lin1_c);

% Diagrama de fase para el sistema linealizado alrededor de Xeq1
figure(4); 

plot3(X1a(:, 1), X1a(:, 2), X1a(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(X1b(:, 1), X1b(:, 2), X1b(:, 3), 'b--', 'LineWidth', 1.5);
plot3(X1c(:, 1), X1c(:, 2), X1c(:, 3), 'r:', 'LineWidth', 1.5);

% Marco el origen, que es el punto xeq1 en este sistema
plot3(0, 0, 0, 'k*', 'MarkerSize', 15, 'LineWidth', 2); 
title('Diagrama de Fase - Sistema Linealizado alrededor de Xeq1');
xlabel('x'); ylabel('y'); zlabel('z');
legend('Cond. iniciales xi_0', 'Cond. iniciales xi_1', 'Cond. iniciales xi_2', 'Punto de Equilibrio Xeq1');
grid on; hold off;
axis equal;

%PUNTO DE EQUILIBRIO 2
% (perturbación = C.I. originales - punto de equilibrio)
X0_lin2_a = xi_0_col - xeq2; % C.I. 1 (desde xi_0)
X0_lin2_b = xi_1_col - xeq2; % C.I. 2 (desde xi_1)
X0_lin2_c = xi_2_col - xeq2; % C.I. 3 (desde xi_2)

% Simulación de las 3 trayectorias para diferentes CI, para el sistema linealizado f_lin2
[t2a, X2a] = ode45(f_lin2, t_span, X0_lin2_a);
[t2b, X2b] = ode45(f_lin2, t_span, X0_lin2_b);
[t2c, X2c] = ode45(f_lin2, t_span, X0_lin2_c);

figure;

% Subplot 1: Trayectoria desde xi_0
subplot(1, 3, 1);
plot3(X2a(:, 1), X2a(:, 2), X2a(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(0, 0, 0, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
title('Desde Cond. Inicial xi_0');
xlabel('x'); ylabel('y'); zlabel('z');
legend('Trayectoria', 'Xeq2');
grid on; axis equal; hold off;

% Subplot 2: Trayectoria desde xi_1
subplot(1, 3, 2);
plot3(X2b(:, 1), X2b(:, 2), X2b(:, 3), 'b--', 'LineWidth', 1.5); hold on;
plot3(0, 0, 0, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
title('Desde Cond. Inicial xi_1');
xlabel('x'); ylabel('y'); zlabel('z');
legend('Trayectoria', 'Xeq2');
grid on; axis equal; hold off;

% Subplot 3: Trayectoria desde xi_2
subplot(1, 3, 3);
plot3(X2c(:, 1), X2c(:, 2), X2c(:, 3), 'r:', 'LineWidth', 1.5); hold on;
plot3(0, 0, 0, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
title('Desde Cond. Inicial xi_2');
xlabel('x'); ylabel('y'); zlabel('z');
legend('Trayectoria', 'Xeq2');
grid on; axis equal; hold off;

%% Soluciones temporales para modelo Linealizado 

% Soluciones Temporales para Xeq1 - Usando Condicion inicial xi_0 
figure(6); 
subplot(3,1,1);
plot(t1a, X1a(:,1), 'g');
xlabel('Tiempo');
ylabel('x(t)');
title('Soluciones Temporales Linealizado en Xeq1 - Condicion Inicial xi_0');
grid on;

subplot(3,1,2);
plot(t1a, X1a(:,2), 'b');
xlabel('Tiempo');
ylabel('y(t)');
grid on;

subplot(3,1,3);
plot(t1a, X1a(:,3), 'r');
xlabel('Tiempo');
ylabel('z(t)');
grid on;

% Soluciones Temporales para Xeq2 - Usando Condicion inicial xi_0
figure(7); 
subplot(3,1,1);
plot(t2a, X2a(:,1), 'g');
xlabel('Tiempo');
ylabel('x(t)');
title('Soluciones Temporales Linealizado en Xeq2 - Condición Inicial xi_0');
grid on;

subplot(3,1,2);
plot(t2a, X2a(:,2), 'b');
xlabel('Tiempo');
ylabel('y(t)');
grid on;

subplot(3,1,3);
plot(t2a, X2a(:,3), 'r');
xlabel('Tiempo');
ylabel('z(t)');
grid on;

%% Función Transferencia

% primer punto de eq
% matrices sistema de estados - A,B,C,D

%comienzo con 1 0 0 para observar la salida de X
A1 = J_eval1; %xeq1
A2 = J_eval2; %xeq2
B = [1; 0; 0]; % B=  matriz de entrada del sistema en espacio de estados -> matriz columnade 3x1 = el sistema tiene 3 estados
C = [1 0 0]; % C = matriz salida 
D = 0; % matriz de transmisión directa del sistema -> como no hay se pone 0

[num1, den1] = ss2tf(A1, B, C, D, 1); %convierte un sistema en espacio de estados a una función de transferencia //resul en vectores
[num2, den2] = ss2tf(A2, B, C, D, 1);


%%%
% disp('XEQ1' \n,'Coeficientes del Numerador:');
%disp(num1);
%disp('Coeficientes del Denominador:');
%disp(den1);


%disp('XEQ2' \n,'Coeficientes del Numerador:');
%disp(num2);
%disp('Coeficientes del Denominador:');
%disp(den2);


G1= tf(num1, den1); %crea la función de transferencia a partir de los coeficientes del numerador y denominador obtenidos.

%disp('Función de Transferencia G1(s) [U(s) -> X(s)] para Xeq1:');
%disp(G1_tf);

G2= tf(num2, den2);

%disp('Función de Transferencia G2(s) [U(s) -> X(s)] para Xeq2:');
%disp(G2_tf);

%Y
B = [0; 1; 0]; 
C = [0 1 0];
[num1, den1] = ss2tf(A1, B, C, D); 
[num2, den2] = ss2tf(A2, B, C, D); 
H1 = tf(num1, den1); % Función de transferencia
H2 = tf(num2, den2); % Función de transferencia

%CALCULO PARA Z
B = [0; 0; 1];
C = [0 0 1]; 
[num1, den1] = ss2tf(A1, B, C, D); 
[num2, den2] = ss2tf(A2, B, C, D); 
L1 = tf(num1, den1);
L2 = tf(num2, den2); 

%% Respuesta al Impulso 

%Graficos
figure(8);
title('Impulso - Linealizado en Xeq1');
subplot(3,1,1);
step(G1, t);
title('Respuesta al Escalón en X');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,2);
step(H1, t);
title('Respuesta al Escalón en Y');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,3);
step(L1, t);
title('Respuesta al Escalón en Z');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

figure(9)
subplot(3,1,1);
step(G2, t);
title('Respuesta al Escalón en X');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,2);
step(H2, t);
title('Respuesta al Escalón en Y');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,3);
step(L2, t);
title('Respuesta al Escalón en Z');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

%% Respuesta al escalón
figure(10)
subplot(3,1,1);
step(G1, t);
title('Respuesta al Escalón X');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,2);
step(H1, t);
title('Respuesta al Escalón Y');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,3);
step(L1, t);
title('Respuesta al Escalón Z');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

%X2
figure(11)
subplot(3,1,1);
step(G2, t);
title('Respuesta al Escalón variable X');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,2);
step(H2, t);
title('Respuesta al Escalón variable Y');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;
subplot(3,1,3);
step(L2, t);
title('Respuesta al Escalón variable Z');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;


%% Diagrama polos y ceros
figure(12)

set(gcf, 'Position', [100, 100, 1200, 400]) % Ajusta tamaño de la figura

% Colores y marcadores
colorPolos = 'r'; % rojo para polos
colorCeros = 'b'; % azul para ceros
markerPolos = 'x';
markerCeros = 'o';

% Subplot para X
subplot(1,3,1)
[pG,zG] = pzmap(G1);
plot(real(pG), imag(pG), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zG), imag(zG), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('X (xeq1)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Subplot para Y
subplot(1,3,2)
[pH,zH] = pzmap(H1);
plot(real(pH), imag(pH), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zH), imag(zH), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('Y', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Subplot para Z
subplot(1,3,3)
[pL,zL] = pzmap(L1);
plot(real(pL), imag(pL), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zL), imag(zL), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('Z', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Título general
sgtitle('Diagramas de Polos y Ceros para xeq1 en X, Y y Z', 'FontSize', 14, 'FontWeight', 'bold');



figure(13)
set(gcf, 'Position', [100, 100, 1200, 400]) % Ajusta tamaño de la figura

% Subplot para X
subplot(1,3,1)
[pG,zG] = pzmap(G2);
plot(real(pG), imag(pG), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zG), imag(zG), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('X (xeq1)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Subplot para Y
subplot(1,3,2)
[pH,zH] = pzmap(H2);
plot(real(pH), imag(pH), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zH), imag(zH), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('Y', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Subplot para Z
subplot(1,3,3)
[pL,zL] = pzmap(L2);
plot(real(pL), imag(pL), markerPolos, 'Color', colorPolos, 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(zL), imag(zL), markerCeros, 'Color', colorCeros, 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
title('Z', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje Real'); ylabel('Eje Imaginario');
legend('Polos', 'Ceros', 'Location', 'best');

% Título general
sgtitle('Diagramas de Polos y Ceros para xeq2 en X, Y y Z', 'FontSize', 14, 'FontWeight', 'bold');
