clc
clear
close all
n=0;
%% Parámetros del sistema
a = 0.2; b = 0.2; c = 5.7;

% Condiciones iniciales
xi_0 = [0.01, -0.03, 0.04]; %antes del primer punto de eq
xi_1 = [6, -28, 28];%desp del 1ro antes del segundo
xi_2 = [1 1 1]; %despues del segundo


% Ventana temporal
t_span = linspace(0, 200, 5000);

% Sistema de Rössler

rossler = @(t,x) [-x(2) - x(3);
           x(1) + a*x(2);
           b + x(3)*(x(1) - c)];

% Resolución
[t, cond1] = ode45(rossler, t_span, xi_0);
[~, cond2] = ode45(rossler, t_span, xi_1);
[~, cond3] = ode45(rossler, t_span, xi_2);

%% TRAYECTORIAS TEMPORALES SISTEMA NO LINEALIZADO
n=n+1;
figure(n);
plot3(cond1(:, 1), cond1(:, 2), cond1(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(cond2(:, 1), cond2(:, 2), cond2(:, 3), 'b--', 'LineWidth', 1.5);
plot3(cond3(:, 1), cond3(:, 2), cond3(:, 3), 'r-', 'LineWidth', 1.5);
plot3(xi_0(1), xi_0(2), xi_0(3), 'r*', 'MarkerSize',10,'MarkerFaceColor', 'g'); % Inicio X01
plot3(xi_1(1), xi_1(2), xi_1(3), 'go', 'MarkerFaceColor', 'b'); % Inicio x1
plot3(xi_2(1), xi_2(2), xi_2(3), 'bo', 'MarkerFaceColor', 'r'); % Inicio x2
xlabel('x'); ylabel('y'); zlabel('z');
title('Diagrama de Fase del Sistema No Lineal');
legend('Trayectoria 1', 'Trayectoria 2', 'Trayectoria 3', 'Inicio 1', 'Inicio 2', 'Inicio 3');
grid on;hold off;

% Gráfico de trayectorias temporales sistema no linealizado PRIMER PUNTO INICIAL
n=n+1;
figure(n);
subplot(3,1,1)
plot(t, cond1(:,1), 'g'); 
plot(t, cond2(:,1), 'r'); 
plot(t, cond3(:,1), 'b'); 
xlabel('Tiempo');
ylabel('x(t)');
xlim([0 200]);
title('Trayectorias Temporales de las soluciones ODE - No Lineal');
grid on
subplot(3,1,2)
plot(t, cond1(:,2), 'b'); 
xlabel('Tiempo');
ylabel('y(t)');
xlim([0 200]);
grid on
subplot(3,1,3)
plot(t, cond1(:,3), 'r'); 
ylabel('z(t)');
xlim([0 200]);
grid on;


% Gráfico de trayectorias temporales sistema no linealizado SEGUNDO PUNTO INICIAL
n=n+1;
figure(n);
subplot(3,1,1)
plot(t, cond2(:,1), 'g'); 
xlabel('Tiempo');
ylabel('x(t)');
xlim([0 200]);
title('Trayectorias Temporales de las soluciones ODE - No Lineal');
grid on
subplot(3,1,2)
plot(t, cond2(:,2), 'b'); 
xlabel('Tiempo');
ylabel('y(t)');
xlim([0 200]);
grid on
subplot(3,1,3)
plot(t, cond2(:,3), 'r'); 
ylabel('z(t)');
xlim([0 200]);
grid on;


% Gráfico de trayectorias temporales sistema no linealizado TERCER PUNTO INICIAL
n=n+1;
figure(n);
subplot(3,1,1)
plot(t, cond3(:,1), 'g'); 
xlabel('Tiempo');
ylabel('x(t)');
xlim([0 200]);
title('Trayectorias Temporales de las soluciones ODE - No Lineal');
grid on
subplot(3,1,2)
plot(t, cond3(:,2), 'b'); 
xlabel('Tiempo');
ylabel('y(t)');
xlim([0 200]);
grid on
subplot(3,1,3)
plot(t, cond3(:,3), 'r'); 
ylabel('z(t)');
xlim([0 200]);
grid on;


%% prube chat
n = n + 1;
figure(n);

% Ajustar tamaño y fuente
set(gcf, 'Position', [100, 100, 800, 600]); % tamaño ventana
sgtitle('Trayectorias Temporales de las soluciones ODE - No Lineal para 3 condiciones iniciales', ...
        'FontSize', 20, 'FontWeight', 'bold'); % título general

% Paleta de colores
col1 = [0 0.6 0];    % Verde intenso
col2 = [0 0 0.8];    % Azul profundo
col3 = [0.8 0 0];    % Rojo fuerte

% Subplot 1: X(t)
subplot(3,1,1);
plot(t, cond1(:,1), 'Color', col1, 'LineWidth', 0.5); hold on;
plot(t, cond2(:,1), 'Color', col2, 'LineWidth', 0.5);
plot(t, cond3(:,1), 'Color', col3, 'LineWidth', 0.5);
xlabel('Tiempo', 'FontSize', 12);
ylabel('x(t)', 'FontSize', 12);
title('Variación de X', 'FontSize', 12);
xlim([0 200]);
legend('x01', 'x02', 'x03', 'Location', 'best');
grid on;

% Subplot 2: Y(t)
subplot(3,1,2);
plot(t, cond1(:,2), 'Color', col1, 'LineWidth', 0.5); hold on;
plot(t, cond2(:,2), 'Color', col2, 'LineWidth', 0.5);
plot(t, cond3(:,2), 'Color', col3, 'LineWidth', 0.5);
xlabel('Tiempo', 'FontSize', 12);
ylabel('y(t)', 'FontSize', 12);
title('Variación de Y', 'FontSize', 12);
xlim([0 200]);
grid on;

% Subplot 3: Z(t)
subplot(3,1,3);
plot(t, cond1(:,3), 'Color', col1, 'LineWidth', 0.5); hold on;
plot(t, cond2(:,3), 'Color', col2, 'LineWidth', 0.5);
plot(t, cond3(:,3), 'Color', col3, 'LineWidth', 0.5);
xlabel('Tiempo', 'FontSize', 12);
ylabel('z(t)', 'FontSize', 12);
title('Variación de Z', 'FontSize', 12);
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
%% prube chat
% Supongamos que ya tienes cond1, cond2, cond3 definidos (trayectorias)
% y xi_0, xi_1, xi_2 como puntos iniciales.
n=n+1;
figure(n);
set(gcf,'Color','w'); % Fondo blanco

% Colores para cada condición
colors = {'g','b','r'}; % verde, azul, rojo

% --- Subplot 1: Plano XY ---
subplot(2,2,1);
plot(cond1(:,1), cond1(:,2), 'Color', colors{1}, 'LineWidth', .5); hold on;
plot(cond2(:,1), cond2(:,2), 'Color', colors{2}, 'LineWidth', .5);
plot(cond3(:,1), cond3(:,2), 'Color', colors{3}, 'LineWidth', .5);
plot(xi_0(1), xi_0(2), '*', 'Color', colors{1}, 'MarkerSize', 8);
plot(xi_1(1), xi_1(2), 'o', 'Color', colors{2}, 'MarkerSize', 6);
plot(xi_2(1), xi_2(2), 'o', 'Color', colors{3}, 'MarkerSize', 6);
legend('Trayectoria 1', 'Trayectoria 2', 'Trayectoria 3', 'Inicio 1', 'Inicio 2', 'Inicio 3','FontSize', 16);
xlabel('X'); ylabel('Y'); title('Plano XY'); grid on; axis tight;

% --- Subplot 2: Plano ZX ---
subplot(2,2,2);
plot(cond1(:,3), cond1(:,1), 'Color', colors{1}, 'LineWidth', .5); hold on;
plot(cond2(:,3), cond2(:,1), 'Color', colors{2}, 'LineWidth', .5);
plot(cond3(:,3), cond3(:,1), 'Color', colors{3}, 'LineWidth', .5);
plot(xi_0(3), xi_0(1), '*', 'Color', colors{1}, 'MarkerSize', 8);
plot(xi_1(3), xi_1(1), 'o', 'Color', colors{2}, 'MarkerSize', 6);
plot(xi_2(3), xi_2(1), 'o', 'Color', colors{3}, 'MarkerSize', 6);
xlabel('Z'); ylabel('X'); title('Plano ZX'); grid on; axis tight;

% --- Subplot 3: Plano YZ ---
subplot(2,2,3);
plot(cond1(:,2), cond1(:,3), 'Color', colors{1}, 'LineWidth', .5); hold on;
plot(cond2(:,2), cond2(:,3), 'Color', colors{2}, 'LineWidth', .5);
plot(cond3(:,2), cond3(:,3), 'Color', colors{3}, 'LineWidth', .5);
plot(xi_0(2), xi_0(3), '*', 'Color', colors{1}, 'MarkerSize', 8);
plot(xi_1(2), xi_1(3), 'o', 'Color', colors{2}, 'MarkerSize', 6);
plot(xi_2(2), xi_2(3), 'o', 'Color', colors{3}, 'MarkerSize', 6);
xlabel('Y'); ylabel('Z'); title('Plano YZ'); grid on; axis tight;

% --- Subplot 4: Vista 3D ---
subplot(2,2,4);
plot3(cond1(:,1), cond1(:,2), cond1(:,3), 'Color', colors{1}, 'LineWidth', .5); hold on;
plot3(cond2(:,1), cond2(:,2), cond2(:,3), 'Color', colors{2}, 'LineWidth', .5);
plot3(cond3(:,1), cond3(:,2), cond3(:,3), 'Color', colors{3}, 'LineWidth', .5);
plot3(xi_0(1), xi_0(2), xi_0(3), '*', 'Color', colors{1}, 'MarkerSize', 8);
plot3(xi_1(1), xi_1(2), xi_1(3), 'o', 'Color', colors{2}, 'MarkerSize', 6);
plot3(xi_2(1), xi_2(2), xi_2(3), 'o', 'Color', colors{3}, 'MarkerSize', 6);
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Vista 3D');



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

%% Linearizacion del sistema sin condiciones iniciales
disp('Jacobiana en el punto de equilibrio x1:')
J_eval1 = double(subs(J, [x y z], xeq1'));
disp(J_eval1)

f_lin1 = @(t, X) J_eval1 * X; %declaracion de la funcion 

disp('Jacobiana en el punto de equilibrio x2:')
J_eval2 = double(subs(J, [x y z], xeq2'));
disp(J_eval2)

f_lin2 = @(t, X) J_eval2 * X;


t_span = linspace(0, 50, 2000); %para sistema lineal necesito tiempo más corto

% Simulación
[t, X1] = ode45(f_lin1, t_span, xeq1); %pasar xeq1 ya lo resta
[~, X2] = ode45(f_lin2, t_span, xeq2);

% Gráfico
n=n+1;
figure(n);
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

%% Autovectores

% Calcular autovectores para J_eval1 y J_eval2
[V1, D1] = eig(J_eval1);  % V1 = autovectores, D1 = autovalores en diagonal
[V2, D2] = eig(J_eval2);

disp('Autovectores para el primer punto de equilibrio:')
disp(V1)

disp('Autovectores para el segundo punto de equilibrio:')
disp(V2)


% Método manual usando null()
autovectores1 = zeros(size(J_eval1)); % inicializar
for i = 1:length(lambda1)
    v = null(J_eval1 - lambda1(i)*eye(size(J_eval1)));
    autovectores1(:,i) = v; % guardar autovector
end

autovectores2 = zeros(size(J_eval2));
for i = 1:length(lambda2)
    v = null(J_eval2 - lambda2(i)*eye(size(J_eval2)));
    autovectores2(:,i) = v;
end

disp('Autovectores para J_eval1:')
disp(autovectores1)

disp('Autovectores para J_eval2:')
disp(autovectores2)



%% Matriz de transicion de estados


%% Matriz fundamental

%% Simulaciones para modelo linealizado - Diagramas de Fase

% Repetición simulaciones con condiciones iniciales originales

% Llevar condiciones iniciales al sistema de coordenadas de xeq1 - mover
% centro de referencia

% condiciones iniciales a partir de mi nuevo cero (perturbación = C.I. originales - punto de equilibrio)
X0_lin1_a = xi_0' - xeq1; % C.I. 1 (desde xi_0)
X0_lin1_b = xi_1' - xeq1; % C.I. 2 (desde xi_1)
X0_lin1_c = xi_2' - xeq1; % C.I. 3 (desde xi_2)

% Simulación de las 3 trayectorias para diferentes CI, para el sistema linealizado f_lin1
[t1a, X1a] = ode45(f_lin1, t_span, X0_lin1_a);
[t1b, X1b] = ode45(f_lin1, t_span, X0_lin1_b);
[t1c, X1c] = ode45(f_lin1, t_span, X0_lin1_c);

% Diagrama de fase para el sistema linealizado alrededor de Xeq1
n=n+1;
figure(n);

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
X0_lin2_a = xi_0' - xeq2; % C.I. 1 (desde xi_0)
X0_lin2_b = xi_1' - xeq2; % C.I. 2 (desde xi_1)
X0_lin2_c = xi_2' - xeq2; % C.I. 3 (desde xi_2)

% Simulación de las 3 trayectorias para diferentes CI, para el sistema linealizado f_lin2
[t2a, X2a] = ode45(f_lin2, t_span, X0_lin2_a);
[t2b, X2b] = ode45(f_lin2, t_span, X0_lin2_b);
[t2c, X2c] = ode45(f_lin2, t_span, X0_lin2_c);

% Diagrama de fase para el sistema linealizado alrededor de Xeq1
n=n+1;
figure(n);

plot3(X2a(:, 1), X2a(:, 2), X2a(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(X2b(:, 1), X2b(:, 2), X2b(:, 3), 'b--', 'LineWidth', 1.5);
plot3(X2c(:, 1), X2c(:, 2), X2c(:, 3), 'r:', 'LineWidth', 1.5);

% Marco el origen, que es el punto xeq1 en este sistema
plot3(0, 0, 0, 'k*', 'MarkerSize', 15, 'LineWidth', 2); 
title('Diagrama de Fase - Sistema Linealizado alrededor de Xeq2');
xlabel('x'); ylabel('y'); zlabel('z');
legend('Cond. iniciales xi_0', 'Cond. iniciales xi_1', 'Cond. iniciales xi_2', 'Punto de Equilibrio Xeq2');
grid on; hold off;
axis equal;

%% Soluciones para sistema linealizado

% Convertir soluciones linealizadas a coordenadas originales
X1a_orig = X1a + xeq1';
X1b_orig = X1b + xeq1';
X1c_orig = X1c + xeq1';

X2a_orig = X2a + xeq2';
X2b_orig = X2b + xeq2';
X2c_orig = X2c + xeq2';

% -------- Xeq1 --------
n = n + 1;
figure(n);

% x(t)
subplot(3,1,1);
plot(t1a, X1a_orig(:,1), 'g'); hold on;
plot(t1a, X1b_orig(:,1), 'r');
plot(t1a, X1c_orig(:,1), 'b');
xlabel('Tiempo'); ylabel('x(t)');
title('Soluciones Temporales Linealizado en Xeq1 (coordenadas originales)');
legend('xi_0','xi_1','xi_2','Location','best'); grid on;
xlim([0 5]);   % <-- aplicar al final

% y(t)
subplot(3,1,2);
plot(t1a, X1a_orig(:,2), 'g'); hold on;
plot(t1a, X1b_orig(:,2), 'r');
plot(t1a, X1c_orig(:,2), 'b');
xlabel('Tiempo'); ylabel('y(t)'); grid on;
legend('xi_0','xi_1','xi_2','Location','best');
xlim([0 5]);   % <-- aplicar al final

% z(t)
subplot(3,1,3);
plot(t1a, X1a_orig(:,3), 'g'); hold on;
plot(t1a, X1b_orig(:,3), 'r');
plot(t1a, X1c_orig(:,3), 'b');
xlabel('Tiempo'); ylabel('z(t)'); grid on;
legend('xi_0','xi_1','xi_2','Location','best');
xlim([0 5]);   % <-- aplicar al final


% -------- Xeq2 --------
n = n + 1;
figure(n);

% x(t)
subplot(3,1,1);
plot(t2a, X2a_orig(:,1), 'g'); hold on;
plot(t2a, X2b_orig(:,1), 'r');
plot(t2a, X2c_orig(:,1), 'b');
xlabel('Tiempo'); ylabel('x(t)');
title('Soluciones Temporales Linealizado en Xeq2 (coordenadas originales)');
legend('xi_0','xi_1','xi_2','Location','best'); grid on;
xlim([0 5]);   % <-- aplicar al final

% y(t)
subplot(3,1,2);
plot(t2a, X2a_orig(:,2), 'g'); hold on;
plot(t2a, X2b_orig(:,2), 'r');
plot(t2a, X2c_orig(:,2), 'b');
xlabel('Tiempo'); ylabel('y(t)'); grid on;
legend('xi_0','xi_1','xi_2','Location','best');
xlim([0 5]);   % <-- aplicar al final

% z(t)
subplot(3,1,3);
plot(t2a, X2a_orig(:,3), 'g'); hold on;
plot(t2a, X2b_orig(:,3), 'r');
plot(t2a, X2c_orig(:,3), 'b');
xlabel('Tiempo'); ylabel('z(t)'); grid on; 
legend('xi_0','xi_1','xi_2','Location','best');
xlim([0 5]);
%% Función Transferencia

%comienzo con 1 0 0 para observar la salida de X
A1 = J_eval1; %xeq1
A2 = J_eval2; %xeq2
B = [0 1; 1 0; 0 0]; % B =  matriz de entrada del sistema en espacio de estados -> matriz columnade 3x1 = el sistema tiene 3 estados
C = [1 0 0]; % C = matriz salida 
D = [0 0]; % matriz de transmisión directa del sistema -> como no hay se pone 0

[num1_in1, den1] = ss2tf(A1, B, C, D, 1); 
[num1_in2, ~]   = ss2tf(A1, B, C, D, 2); % El 2 es porque es para la entrada 2, como B es 3x2 tengo 2 entradas

[num2_in1, den2] = ss2tf(A2, B, C, D, 1); 
[num2_in2, ~]   = ss2tf(A2, B, C, D, 2); 

% Crea la función de transferencia a partir de los coeficientes del numerador y denominador obtenidos.
% x1eq
G1_in1 = tf(num1_in1, den1); % Entrada 1
G1_in2 = tf(num1_in2, den1); % Entrada 2
%disp('Función de Transferencia G1(s) [U(s) -> X(s)] para Xeq1:');
%disp(G1_tf);

% x2eq
G2_in1 = tf(num2_in1, den2);
G2_in2 = tf(num2_in2, den2);

%disp('Función de Transferencia G2(s) [U(s) -> X(s)] para Xeq2:');
%disp(G2_tf);

%% Respuesta al Impulso 

%Graficos
n=n+1;
figure(n);
subplot(2,2,1);
impulse(G1_in1, t);
title('Respuesta al Impulso - Entrada 1');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

subplot(2,2,2);
impulse(G1_in2, t);
title('Respuesta al Impulso - Entrada 2');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

sgtitle('Sistema linealizado alrededor de Xeq1')

n=n+1;
figure(n);
subplot(2,2,1);
impulse(G2_in1, t);
title('Respuesta al Impulso - Entrada 1');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

subplot(2,2,2);
impulse(G2_in2, t);
title('Respuesta al Impulso - Entrada 2');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

sgtitle('Sistema linealizado alrededor de Xeq2')


%% Respuesta al escalon

%Graficos
n=n+1;
figure(n);
subplot(2,2,1);
step(G1_in1, t);
title('Respuesta al Escalon - Entrada 1');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

subplot(2,2,2);
step(G1_in2, t);
title('Respuesta al Escalon - Entrada 2');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

sgtitle('Sistema linealizado alrededor de Xeq1')

n=n+1;
figure(n);

subplot(2,2,1);
step(G2_in1, t);
title('Respuesta al Escalon - Entrada 1');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

subplot(2,2,2);
step(G2_in2, t);
title('Respuesta al Escalon - Entrada 2');
xlabel('Tiempo');
ylabel('Amplitud');
grid on;

sgtitle('Sistema linealizado alrededor de Xeq2')