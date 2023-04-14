% trip.m
% for AE68714 (Pusan National Univeristy)
%
% Periodic tridiagonal solver using Sherman-Morrison algorithm
% --                -- --------
% | b_1 c_1              a_1 | | fout_1 | = | fin_1 |
% | a_2 b_2 c_2              | | fout_2 | = | fin_2 |
% |     a_3 b_3 c_3          | | fout_3 | = | fin_3 |
% |      .    .   .          | |   .    | = |   .   |
% |           .   .   .      | |   .    | = |   .   |
% |               .   .   .  | |   .    | = |   .   |
% | c_M              a_M b_M | | fout_M | = | fin_M |
% ---- --------
% Note: x and fdum are dummy arrays
%f array is overwritten by solution
%a, b and c are preserved
function f = trip(a,b,c,f)
M = length(a);
x = zeros(M,1);
% modified forcing term for second tridiagonal
fdum = zeros(M,1);
fdum(1) = -a(1);
fdum(M-1) = -c(M-1);

% solve two tridiagonals
x(1) = c(1)/b(1);
f(1) = f(1)/b(1);
fdum(1) = fdum(1)/b(1);

% forward sweep for two systems
% Note: can simplify fdum since all zero initially from ji=2:m-2
for ji=2:M-2
z = 1./( b(ji) - a(ji)*x(ji-1) );
x(ji) = c(ji)*z;
f(ji) = ( f(ji) - a(ji)*f(ji-1) )*z;
fdum(ji) = - a(ji)*fdum(ji-1)*z;   % last culume a1 move from up to down
end

% need to include since fdum(M-1) not initially zero
ji = M-1;
z = 1./( b(ji) - a(ji)*x(ji-1) );
x(ji) = c(ji)*z;
f(ji) = ( f(ji) - a(ji)*f(ji-1) )*z;
fdum(ji) = ( fdum(ji) - a(ji)*fdum(ji-1) )*z; 

% backward sweep for two systems
for ji=M-2:-1:1
f(ji) = f(ji)- x(ji)*f(ji+1);
fdum(ji) = fdum(ji)- x(ji)*fdum(ji+1);
end
% now recombine the two solutions
% solve at the last point
z = 1./ ( b(M)+a(M)*fdum(M-1)+c(M)*fdum(1) );
f(M) = ( f(M)-a(M)*f(M-1)-c(M)*f(1) )*z;
% now solve the rest of the points
for ji=1:M-1
f(ji) = f(ji) + fdum(ji)*f(M);
end
end