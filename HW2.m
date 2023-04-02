syms k 
x = linspace(0, 2*pi); %setup x

%{
u = symsum((1/k)*sin(2*k*x), k, 1, 6); %u(x)
u1 = symsum((1/k)*sin(2*k*x), k, 1, 1000); %u(x) inf terms
%}


u2 = symsum(2*cos(2*k*x), k, 1, 6); %u'(x) 6 terms
u3 = symsum(2*cos(2*k*x), k, 1, 1000); %u'(x) inf terms

%plot u(x) and u1(x)
%{
plot(x, u, 'LineWidth',2);
hold on
plot(x, u1, 'LineWidth',2);
hold on
%}

%plot u'(x) and u'1(x)
plot(x, u2, 'LineWidth',2);
hold on
plot(x, u3, 'LineWidth',2);

xlabel('x')
ylabel('u''(x)')
%legend('6 Terms', 'inf terms')
legend('1st derivative (6 terms)', '1st derivative (inf terms)')

%{
syms x
truedif = diff(symsum((1/k)*sin(2*k*x), k, 1, 6)); %1st derivative of u(x)

%forward dt=0.04
uj_foward1 = symsum((1/k) *sin(2*k*(pi/8 + 0.04)), k, 1, 6); 
%backward dt=0.04
uj_backward1 = symsum((1/k) *sin(2*k*(pi/8 - 0.04)), k, 1, 6);

%forward dt=0.02
uj_foward2 = symsum((1/k) *sin(2*k*(pi/8 + 0.02)), k, 1, 6);
%backward dt=0.02
uj_backward2 = symsum((1/k) *sin(2*k*(pi/8 - 0.02)), k, 1, 6);

%approximation 1st derivative dt=0.04 three points central
duj1 = (uj_foward1-uj_backward1)/(2*0.04);

%approximation 1st derivative dt=0.02 three points central
duj2 = (uj_foward2-uj_backward2)/(2*0.02);

err1 = -3.41421 - duj1; %Error for dt=0.04
err2 = -3.41421 - duj2; %Error for dt=0.02
%}