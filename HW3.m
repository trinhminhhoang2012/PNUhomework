% Code to demonstrate derivative accuracy.
% for AE68714 (Pusan National Univeristy)
%
% Using a periodic made up function 
% u = sum 1./k*sin(k*(x)); k = 1, 6
%
% Compares its first derivative with numerical derivatives computed
% using several difference methods:
% 1) 2nd order Central du/dx = B(-1,0,1)/(2*dx) u
% 2) 3rd order Backwards Biased du/dx = B(1,-6,3,2,0)/(6*dx) u
%
% Code runs through a sequence of mesh refinements, halving dx with each.
% Error is defined as the L1 norm of the difference between 
% the exact derivative and numerical difference.
%
clear;
for ifig=1:7
    figure(ifig);
clf;

end
% Let's set some things here and ask for input
% CHANGE THIS IF ADD MORE METHODS
num_methods = 2; % total number of approx. methods to look at
% Ask for input
num_levels = input('\nNumber of levels (<=4) of grid refinement = ? ');
method_plot = input(['Which method to plot? (<=',num2str(num_methods),') = ? ']);
% Set these up for the various levels of mesh refinement and methods
DX = zeros(num_levels,1);
ERROR = zeros(num_levels,num_methods);
%
% Now let's loop through the mesh size levels
%
for level = 1:num_levels
    % This sets jmax = 16 for level=1; increases by factor of two at each level
    jmax = 16*2^(level-1); % >>> adjust the number of mesh point? 
    fprintf('\n For this level: jmax = %5.0f \n',jmax);
    
    % Initialize array sizes for this level
    jm= zeros(jmax,1);
    jp= zeros(jmax,1);
    u = zeros(jmax,1);
    du_exact = zeros(jmax,1);
    du_approx = zeros(jmax,num_methods);
    error = zeros(jmax,num_methods);
    
    % Grid spacing and mesh for this level; note x starts at dx (not zero!)
    dx = 2.*pi/(jmax);
    x = dx*(1:jmax)';
    DX(level) = dx;
    
    % Arrays for periodicity on this level
    for j=1:jmax
        jp(j) = j+1;
        jm(j) = j-1;
    end
    jp(jmax)=1;
    jm(1)=jmax;
    
    % Exact function and derivative on this level
    kbeg=1; kend=6;
 
    for k=kbeg:kend
        u = u + 1./k*sin(2.*k*x);
        du_exact = du_exact + 2.*cos(2.*k*x);
    end
   
    % Individual point to examine (pi/8 = 2*pi/16) and print out
    j_point = jmax/16;
    fprintf(' du_exact(%12.4e) = %14.8e \n',x(j_point),du_exact(j_point));
    
    
    % Let's calculate the derivatives (using vector notation)
    % CHANGE HERE IF ADD MORE APPROXIMATIONS!
    %
    % 1st method) 2nd order Central
    du_approx(1:jmax,1) = (u(jp)-u(jm))/(2*dx);
    
    % 2nd method) 3rd order Backward
    du_approx(1:jmax,2) = (2*u(jp)+3*u-6*u(jm)+u(jm(jm)))/(6*dx);

    % CHANGE THIS IF ADD MORE METHODS
    % 3rd method) 
    % Note: for Compact schemes, a periodic tridiagonal solver used where
    % solve like a Lagrangian predictor with Hermetian corrector
 

    % Compute local error for each approximation as well as L1 Norm
    for nn=1:num_methods
        error(1:jmax,nn) = du_exact(1:jmax) - du_approx(1:jmax,nn);
        ERROR(level,nn) = norm(error(1:jmax,nn)/jmax,1);
    end
    
    % Print out local approximation and error
    fprintf(' du_approx(%12.4e) = %15.8e \n',x(j_point),du_approx(j_point,method_plot));
    fprintf(' error(pi/8) = %14.8e \n',error(j_point,method_plot));
    
    
    % Now let's make some plots
    %
    % Plot function on this level on first
    figure(1);
    plot(x,u,'r-',x,u,'ro','LineWidth',2.0);
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['Function on mesh with ',num2str(jmax),' points']);
    xlabel('X');
    ylabel('U');
   
    % Plot derivative (exact) on this level on second plot 
    figure(2);
    plot(x,du_exact,'r-',x,du_exact,'ro','LineWidth',2.0);
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['First Derivative of Function on mesh with ',num2str(jmax),' points']);
    xlabel('X');

    ylabel('U_x');
    
    % Plot local errors (absolute) for all methods at this level on third plot
    figure(3);
    semilogy(x,abs(error),'o');
    hold on;
    semilogy(x,abs(error),'-','LineWidth',2.0);
    hold off;
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['Absolute Errors on mesh with ',num2str(jmax),' points']);
    xlabel('X');
    ylabel('ABS(Error)');
    
    % Plot local error (absolute) for particular method on fourth plot
    figure(4);
    if(level==1) 
        semilogy(x,abs(error(:,method_plot)),'m',x,abs(error(:,method_plot)),'mo','LineWidth',2.0); 
    end
 
    if(level==2) 
         semilogy(x,abs(error(:,method_plot)),'r',x,abs(error(:,method_plot)),'rx','LineWidth',2.0); 
    end
    
    if(level==3) 
        semilogy(x,abs(error(:,method_plot)),'g',x,abs(error(:,method_plot)),'g+','LineWidth',2.0); 
    end
 
    if(level==4) 
        semilogy(x,abs(error(:,method_plot)),'b','LineWidth',1.0); 
    end
 
    hold on;
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['Absolute Errors for method ',num2str(method_plot)]);
    xlabel('X');
    ylabel('Error');
    
    % Plot local error for particular method on fifth plot
    figure(5);
    hold on;

    if(level==1) 
        plot(x,error(:,method_plot),'m',x,error(:,method_plot),'mo','LineWidth',2.0); 
    end
 
    if(level==2) 
        plot(x,error(:,method_plot),'r',x,error(:,method_plot),'rx','LineWidth',2.0); 
    end
 
    if(level==3) 
        plot(x,error(:,method_plot),'g',x,error(:,method_plot),'g+','LineWidth',2.0); 
    end
 
    if(level==4) 
        plot(x,error(:,method_plot),'b','LineWidth',1.0); 
    end
 
    hold off;
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['Error for method ',num2str(method_plot)]);
    xlabel('X');
    ylabel('Error');
    
    % Plot local derivative for particular method on sixth plot
    figure(6);
    hold on;

    if(level==1) 
        plot(x,du_approx(:,method_plot),'m',x,du_approx(:,method_plot),'mo','LineWidth',2.0); 
    end
 
    if(level==2) 
        plot(x,du_approx(:,method_plot),'r',x,du_approx(:,method_plot),'rx','LineWidth',2.0); 
    end
 
    if(level==3) 
        plot(x,du_approx(:,method_plot),'g',x,du_approx(:,method_plot),'g+','LineWidth',2.0); 
    end
    
    if(level==4) 
        plot(x,du_approx(:,method_plot),'b','LineWidth',1.0); 
    end
 
    hold off;
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    title(['Derivative Approximation for method ',num2str(method_plot)]);
    xlabel('X');
    ylabel('Approximate U_x');
    
    % Let's pause to digest the results on this mesh level
    fprintf(' HIT ANY KEY TO CONTINUE \n');
    pause;
end

% Now let's plot the global error for all methods for all levels
%
% Plot total Error in log-log. Slope of curves is error order measure.
figure(7);
loglog(DX,ERROR,'o');
hold on;
loglog(DX,ERROR,'-','LineWidth',3.0);
DXT=DX(1:num_levels);
if(num_levels>3) loglog(DXT,[2*2.5^2*DXT.^2,2*2.5^3*DXT.^3,2.5^4*DXT.^4,2.5^6*DXT.^6],'k','LineWidth',1.0); end
hold off;
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(' Error (L1 norm) versus DX ');
xlabel('DX');
ylabel('Error');
fprintf('Done \n \n');
% all finished!
