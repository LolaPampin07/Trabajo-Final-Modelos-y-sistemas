clc
clear all
close all

%% Parámetros

%{
valores iniciales de las variables de estado
intervalo de tiempo considerado
funciones MATLAB utilizadas
gráficos temporales
diagramas de fase en 3D
%}


% Coeficientes tipicos
a = 0.2;
b = 0.2;
c = 5.7;

% Variables de estado
syms x y z

% Condiciones iniciales

% Ventana de tiempo


% Ecuaciones del sistema
f1 = -y - z;
f2 = x + a*y;
f3 = b + z*(x - c);

% Vector de funciones
F = [f1; f2; f3];

%% Puntos de equilibrio (F = 0)
sol = solve(F, [x, y, z]);

x_eq = double(sol.x);
y_eq = double(sol.y);
z_eq = double(sol.z);

xeq1 = [double(sol.x(1)),double(sol.y(1)),double(sol.z(1))];
xeq2 = [double(sol.x(2)),double(sol.y(2)),double(sol.z(2))];

disp('Puntos de equilibrio:')
disp(xeq1)
disp(xeq2)

%% Matriz Jacobiana
J = jacobian(F, [x, y, z]);
disp('Matriz Jacobiana simbólica:')
disp(J)

% Linearizacion del sistema
J_eval1 = subs(J, [x, y, z], xeq1);

% Convertir a numérico
J_eval1_num = double(J_eval1);

disp('Jacobiana en el punto de equilibrio x1:')
for i = 1:size(J_eval1_num,1)
    fprintf('[');
    fprintf(' %.2f', J_eval1_num(i,:)); % 2 decimales
    fprintf(' ]\n');
end


J_eval2 = subs(J, [x, y, z], xeq2);
% Convertir a numérico
J_eval2_num = double(J_eval2);

disp('Jacobiana en el punto de equilibrio x2:')
for i = 1:size(J_eval2_num,1)
    fprintf('[');
    fprintf(' %.2f', J_eval2_num(i,:)); % 2 decimales
    fprintf(' ]\n');
end
  
 %% 4. Autovalores
lambda1 = eig(J_eval1);
lambda2 = eig(J_eval2);
fprintf('Autovalores x1:\n');
for i = 1:length(lambda1)
    fprintf(' %.2f + %.2fi\n', real(lambda1(i)), imag(lambda1(i)));
end
fprintf('Autovalores x2:\n');
for i = 1:length(lambda2)
    fprintf(' %.2f + %.2fi\n', real(lambda2(i)), imag(lambda2(i)));
end

%% Diagrama de fase
%% Grafico de trayectorias
%% Funcion transferencia
