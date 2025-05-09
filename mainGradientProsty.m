clear all 
clc
N = 5;
u = [1 4 2 2 2]; % pierwszy przybliżony ciąg sterowań
K = 15;          % maksymalna liczba iteracji 
eps = 0.2;       % dokładność do kryterium stopu 
c = 0;
x0 = 8;            % stan początkowy

x = zeros(1,N);    % deklaracja wektora stanów
x(1) = x0;
u0 = u;
p = zeros(1,N+1);
b = zeros(1,N);

f = @(x,u) x + 0.5*u;            % równanie stanu

L = @(x,u) (10*(x-3).^2 + u.^2); % funkcja kosztu

J = @(x,u) sum(L(x, u));         % wskaźnik jakości 

Jk = zeros(1,K);
xk = zeros(K,5);
uk = xk;
pn = @(x,p) 20*x -60 +p;
bn = @(u,p) 2*u + 0.5*p;
normBn = norm(bn(u,p(2:N+1)));
k = 1; % numer iteracji
while k <= K & eps < normBn
    for i = 1:(N-1)
        x(i+1) = f(x(i),u(i));
    end
Jk(k) = J(x,u); 
    for i = N:-1:1
        p(i) = pn(x(i), p(i+1));
    end
    for i = 1:N
        b(i) = bn(u(i), p(i+1));
    end
normBn = norm(bn(u,p(2:N+1)));
 
 % Funkcja do minimalizacji J(t)
    J_t = @(t) J_t_function(t, x, u, b, f, L, N);
    
    % Znalezienie optymalnego t za pomocą fminsearch
    options = optimset('Display', 'off'); % Wyłączenie wyświetlania komunikatów
    t_opt = fminsearch(J_t, 0, options); % Szukaj t w okolicy 0
    
    
    % Aktualizacja sterowań u
    u = u - t_opt * b;
    uk(k,:) = u;
    xk(k,:) = x;
    
k = k+1;
end

figure
plot(1:K, Jk, 'o', 'LineStyle','--')
xlabel('Numer iteracji')
ylabel('Wartość wskaźnika wartości')

figure 
plot(0:N-1,xk)
xlabel('Numer stanu n')
ylabel('Wartość stanu x')

figure
plot(0:N-1,u0)
hold on
plot(0:N-1,uk)
xlabel('Numer sterowania n')
ylabel('Wartość sterowania u')
hold off
% Funkcja obliczająca J(t) dla danego t
function cost = J_t_function(t, x, u, b, f, L, N)
    u_new = u - t * b;
    x_new = x;
    for i = 1:(N-1)
        x_new(i+1) = f(x_new(i), u_new(i));
    end
    cost = sum(L(x_new, u_new));
end