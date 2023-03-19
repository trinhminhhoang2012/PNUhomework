%x = [0 pi/3 pi/2 2*pi/3 pi] ;
x = linspace(0, 2*pi);
time = linspace(0, 100);
syms m 
u = symsum((1/m)*sin(m*x), m, 2, 5);

%plot(x, u);
%set(gca,'XLim',[0, pi], 'YLim',[-0.5 0.5]);

ut0 = f(0);
ut1 = f(1);
ut2 = f(2);
ut3 = f(3);
ut4 = f(4);
ut5 = f(5);
plot(x, u);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
hold on
plot(x, ut1);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
hold on
plot(x, ut2);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
hold on
plot(x, ut3);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
hold on
plot(x, ut4);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
hold on
plot(x, ut5);
set(gca, 'XLim', [0 pi], 'YLim', [-0.5 1.5])
legend('u(x,t=0)','u(x,t=1)','u(x,t=2)','u(x,t=3)','u(x,t=4)','u(x,t=5)')
function u3 = f(t)
    syms n k
    v = .04;
    x = linspace(0, pi);
    u1 = symsum((1/n)*exp(-v*(n^2)*t)*sin(n*x), n, 2, 5);
    u3 = u1;
end