clear all 
clc

N = 5;
u = [1 4 2 2 2]; % pierwszy przybliżony ciąg sterowań
K = 15;          % maksymalna liczba iteracji 
eps = 0.2;       % dokładność do kryterium stopu 
x0 = 8;          % stan początkowy

x = zeros(1,N);  % deklaracja wektora stanów
x(1) = x0;
u0 = u;
p = zeros(1,N+1);
b = zeros(1,N);
d = zeros(1,N);  % wektor kierunków sprzężonych

f = @(x,u) x + 0.5*u;            % równanie stanu
L = @(x,u) (10*(x-3).^2 + u.^2); % funkcja kosztu
J = @(x,u) sum(L(x, u));         % wskaźnik jakości 

Jk = zeros(1,K);
xk = zeros(K, N);
uk = zeros(K, N);
pn = @(x,p) 20*x - 60 + p;
bn = @(u,p) 2*u + 0.5*p;

normBn = norm(bn(u, p(2:N+1)));
k = 1; % numer iteracji
b_prev = zeros(1,N); % przechowuje gradient z poprzedniej iteracji

while k <= K && eps < normBn
    % Obliczenie stanów x dla aktualnych sterowań u
    for i = 1:(N-1)
        x(i+1) = f(x(i), u(i));
    end
    Jk(k) = J(x, u); 
    xk(k,:) = x;
    
    % Obliczenie p (sprzężeń)
    for i = N:-1:1
        p(i) = pn(x(i), p(i+1));
    end
    
    % Obliczenie b (gradientu)
    for i = 1:N
        b(i) = bn(u(i), p(i+1));
    end
    normBn = norm(b);
    
    % Wyznaczenie kierunku sprzężonego d
    if k == 1
        d = -b; % W pierwszej iteracji kierunek to -b
    else
        c = (norm(b)^2) / (norm(b_prev)^2); % Współczynnik Fletchera-Reevesa
        d = -b + c * d; % Aktualizacja kierunku sprzężonego
    end
    
    % Funkcja do minimalizacji J(t)
    J_t = @(t) J_t_function(t, x, u, d, f, L, N);
    
    % Znalezienie optymalnego t za pomocą fminsearch
    options = optimset('Display', 'off');
    t_opt = fminsearch(J_t, 0, options);
    
    % Aktualizacja sterowań u w kierunku d
    u = u + t_opt * d; % Uwaga: kierunek d już uwzględnia znak (-b)
    uk(k,:) = u;
    
    % Zapamiętanie gradientu dla następnej iteracji
    b_prev = b;
    
    k = k + 1;
end

% Wykresy
figure
plot(1:k-1, Jk(1:k-1), 'o', 'LineStyle', '--')
xlabel('Numer iteracji')
ylabel('Wartość wskaźnika jakości')

figure 
plot(0:N-1, xk(1:k-1,:))
xlabel('Numer stanu n')
ylabel('Wartość stanu x')

figure
plot(0:N-1, u0)
hold on
plot(0:N-1, uk(1:k-1,:))
xlabel('Numer sterowania n')
ylabel('Wartość sterowania u')
legend show
hold off

% Funkcja obliczająca J(t) dla danego t
function cost = J_t_function(t, x, u, d, f, L, N)
    u_new = u + t * d; % Kierunek d już uwzględnia znak (-b)
    x_new = x;
    for i = 1:(N-1)
        x_new(i+1) = f(x_new(i), u_new(i));
    end
    cost = sum(L(x_new, u_new));
end