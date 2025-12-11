%% TURBOPROP DESING
%% LAB 4: PROPELLER DESING WITH BEMT

%% LAB GROUP J
%% MARIA SAN EMETERIO GONZALEZ  100471680
%% ANA BELEN ROMERA BAENA       100471676
%% WOOJIN LEE                   100576121
%% AIDAN JONES                  100578091

clear all
close all
clc

%% LATEX FORMAT GRAPHS
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% VARIABLES 
% Fixed geometrical parameters
d = 0.25;              % Diameter          [m]
R = d/2;               % radius            [m]
Dh = 0.025;            % Hub diameter      [m]
Rh = Dh/2;             % Hub radius        [m]
m = 4;                 % Number of blades
rho = 1.225;           % Air density
mu = 1.81e-5;          % dynamic viscosity [kg/(ms)]

% Desing conditions
Cl_init = 0.5;         % Cte value of cl
V = 5;                 % Freestream        [m/s]
rpm = 3000;   
n=rpm/60;
w = rpm * 2*pi/60;     % rad/s

% Requirements
T_required = 1.0;      % Required thrust   [N]

%% MORE VARIABLES  
N = 30;
r = linspace(Rh, R, N);
x = r/R;               % nondimensional radius

%% ==================================================================================
%%  FIND VALUE OF V0
% Initial guess from Momentum Theory
% V0_init=v
v0_ini =(-V+sqrt(V^2+4*(T_required*2)/(pi*d^2*rho)))/2; % v was solved using the cuadratic formula for solving 2 degree equations

% find T(V0)-T_req= 0
v0=fzero(@(v0_guess) getT(d,N,n,V,v0_guess,w,Cl_init,r,m,x,rho)-T_required, v0_ini );

V0=v0;

%% =============================================================
%% CALCULATE WITH OPTIMIZED VALUE OF v0 CHORD AND BETA 
for i=1:N
    
    % STEP 1
    gamma_prime(i) = atand((V + V0)/(w * r(i)));
    a(i) = (V0/V) * (cosd(gamma_prime(i)))^2;
    b(i) = (V0/(w * r(i))) * cosd(gamma_prime(i)) * sind(gamma_prime(i));
    
    % STEP 2
    F(i) = 2/pi * acosd(exp(((m * w * d)/(4 * V)) * (x(i) - 1)));
    cir(i) = (pi *x(i)^2 * b(i) * F(i) * w * d^2)/m;
    
    % STEP 3
    alpha_prime = (Cl_init/(2 * pi))*360/(2*pi);
    beta(i) = alpha_prime + gamma_prime(i);
    
    % STEP 4
    Ve(i) = sqrt( (V*(1+a(i)))^2 + ( w*r(i)*(1-b(i)) )^2 );
    s(i) = (2/pi) * ((m * cir(i))/(Cl_init * Ve(i) * x(i) * d));
    % chord calculated to plot it against r
    c(i) = s(i) * (2*pi*r(i)/m);    
end


%% BLADE GEOMETRY PLOT
% twist angle vs r plot 
figure
plot(r,beta,'-b', 'LineWidth', 1.5)
xlabel('$r$')
ylabel('$\beta$')
grid on
title('$\beta~vs~r$')
set(gca, 'FontSize', 12)
saveas(figure(1), 'beta.png')

% Chord vs r plot
figure
plot(r,c,'-b', 'LineWidth', 1.5)
xlabel('$r$')
ylabel('$chord$')
grid on
title('$chord~vs~r$')
set(gca, 'FontSize', 12)
saveas(figure(2), 'chord.png')

%% 3D PLOT + CHORD PLOT WITH THE NACA0012
% import data from the NACA0012 AIRFOIL (x,y)
data=readmatrix('NACAPLOT.txt');
x_air = data(:,1);
y_air = data(:,2);


figure
hold on
% dimension and rotate the NACA0012 airfoil for our chord and twist angle
for i=1:length(c)
    R_matrix = [cosd(beta(i)), -sind(beta(i)); sind(beta(i)), cosd(beta(i))];
    Air_coord = [x_air*c(i)-c(i)/2, y_air*c(i)];
    Air_rot = Air_coord*R_matrix;
    x_v(i,:) = Air_rot(:,1);
    y_v(i,:) = Air_rot(:,2);
    k(i, :) = linspace(r(i), r(i), length(Air_rot(:,1)));
    plot(x_v(i,:), y_v(i,:)) % plot the blade in 2D
end
saveas(figure(3), 'cross_section.png')

% plot the blade in 3D using the previous calculations
figure
surf(x_v, y_v, k,'EdgeColor','none', 'FaceLighting','gouraud');
xlabel('X (Chord)');
ylabel('Y');
zlabel('Z (Radial)');
title('3D Representation of propeller blade');
grid on
axis equal
view(3); 
colormap jet; 
saveas(figure(4), '3d.png')


%% BEMT ==========================================================
% now using the same code as in lab 3 

%% NACA0012
% NACA0012 GIVEN BY PROFESSOR 
% extract the values from the txt file 
CD_NACA0012 = readmatrix('CDNACA0012.txt');
CL_NACA0012 = readmatrix('CLNACA0012.txt');

% locate the values in diferent matrices 
alpha_NACA = CD_NACA0012(:,1);
Cd_NACA=CD_NACA0012(:,2:end);
Cl_NACA=CL_NACA0012(:,2:end);

% Re number from the NACA0012  data file 
Re_NACA0012 = [160000, 360000, 700000, 1000000, 2000000, 5000000];

max_iter = 1e8; % max nummber of iterations

num_J=30;
J = linspace(0.1,1,num_J);

% Inicialize matrices for the results
a_f = zeros(N, num_J);
b_f = zeros(N, num_J);
gamma_prime_f = zeros(N, num_J);
alpha_prime_f = zeros(N, num_J);
Re_f = zeros(N, num_J);
Cl_f = zeros(N, num_J);
Cd_f = zeros(N, num_J);
Ve_f = zeros(N, num_J);

%% BEMT 
for i=1:num_J
    
    % initial values for gamma a and b
    gamma_p = 53;
    a2 = 0.1;
    b2 = 0.01;
    
    % tolerance
    tol = 1e-8;
    
    % start loop 
    for j=1:N
        iter = 0;
        converged = false; % inizialice converged to false so it enter in the while loop
        
        while ~converged && iter <= max_iter
            iter = iter + 1;
            
            % STEP 1
            alpha = beta(j)- gamma_p;
            
            Ve2 = sqrt( ( J(i) * n * d )^2 * (1+a2)^2 + ( w * x(j) * R )^2 * (1-b2)^2 ); % effective velocity
            Re = (rho * Ve2 * c(j) * R) / mu; % Re number
            
            % STEP 2
            [Cl, Cd] = GETClCd(alpha, Re, Re_NACA0012, alpha_NACA, Cl_NACA, Cd_NACA);  % calling the function for calculating cl and d
            
            % STEP 3
            % calculate a and b for next iteration
            a_i = (s(j) / (4 * sind(gamma_p).^2) * (Cl * cosd(gamma_p) - Cd * sind(gamma_p)));
            % a_i=a_cal/1+a_cal
            b_i = (s(j)/ (4 * sind(gamma_p) * cosd(gamma_p)) * (Cl * sind(gamma_p) + Cd * cosd(gamma_p)));
            % b_i=b_cal/1-b_cal
            
            a_cal = a_i/(1-a_i);
            b_cal = b_i/(b_i+1);
            
            % STEP 4
            J_cal = (pi*x(j))*((1-b2)/(1+a2))*tand(gamma_p);
            
            % STEP 5
            if abs(J_cal-J(i))<tol
                converged = true;
            else
                % Change the values for a, b and gamma so it enter the
                % loop again until J_obj is reach
                a2 = a_cal;
                b2 = b_cal;
                gamma_p = gamma_p + 0.1 * (J(i)-J_cal);
            end
        end
        if converged
            % store the calculated values in matrices
            a_f(j,i) = a2;
            b_f(j,i) = b2;
            gamma_prime_f(j,i) = gamma_p;
            alpha_prime_f(j,i) = alpha;
            Re_f(j,i) = Re;
            Cl_f(j,i) = Cl;
            Cd_f(j,i) = Cd;
            Ve_f(j,i) = Ve2;
        else
            % NO CONVERGED
            fprintf('did not converged'); % display message that the iteration did not converged
        end
        
        % STEP 6
        % dCT/dx and dCP/dx
        dC_T(j,i) = (pi^3/4) * s(j) * (Cl_f(j,i) * cosd(gamma_prime_f(j,i)) - Cd_f(j,i) * sind(gamma_prime_f(j,i)))*...
            ((1-b_f(j,i))^2/cosd(gamma_prime_f(j,i))^2) * x(j)^3;
        
        dC_P (j,i)= ((pi^4)/4) * s(j) * (Cl_f(j,i) * sind(gamma_prime_f(j,i)) + Cd_f(j,i) * cosd(gamma_prime_f(j,i)))*...
            ((1-b_f(j,i))^2/cosd(gamma_prime_f(j,i))^2) * x(j)^4;
    end
    % STEP 8
    % Integrate over x using the trapezoidal function from MATLAB
    C_T(i)= trapz(x,dC_T(:,i));
    C_P(i)=trapz(x,dC_P(:,i));
    eta(i)=(C_T(i)/C_P(i)) * J(i);
end

%% BEMT PLOTS
% THRUST COEFFICIENT
figure
plot(J,C_T,'-b', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$C_T$')
grid on
title('$C_T~vs~J$')
set(gca, 'FontSize', 12)
saveas(figure(5), 'ct.png')

% POWER COEFFICIENT
figure
plot(J,C_P,'-r', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$C_P$')
grid on
title('$C_P~vs~J$')
set(gca, 'FontSize', 12)
saveas(figure(6), 'cp.png');

% ETA 
figure
plot(J,eta,'-g', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$\eta$')
grid on
title('$\eta~vs~J$')
set(gca, 'FontSize', 12)
saveas(figure(7), 'eta.png');


%% ============================================
%% FUNCTIONS
% function to obtain the value of the thrust for the optimal value of Vo
% Optimization algorithm andding the alculation of T 
function T = getT(d,N,n,V,v0_guess,w,Cl,r,m,x,rho)
for i=1:N
    
    gamma_prime(i) = atand((V + v0_guess)/(w * r(i)));
    a(i) = (v0_guess/V) * (cosd(gamma_prime(i)))^2;
    b(i) = (v0_guess/(w * r(i))) * cosd(gamma_prime(i)) * sind(gamma_prime(i));
    
    F(i) = 2/pi * acosd(exp(((m * w * d)/(4 * V)) * (x(i) - 1)));
    cir(i) = (pi *x(i)^2 * b(i) * F(i) * w * d^2)/m;
    
    alpha_prime = (Cl/(2 * pi))*360/(2*pi);
    beta(i) = alpha_prime + gamma_prime(i);
    
    Ve(i) = sqrt( (V*(1+a(i)))^2 + ( w*r(i)*(1-b(i)) )^2 );
    s(i) = (2/pi) * ((m * cir(i))/(Cl * Ve(i) * x(i) * d));
    
    c(i) = s(i) * (2*pi*r(i)/m);
    
    Ct_dx(i) = (pi^3/4) * s(i) * Cl*cosd(gamma_prime(i)) * ((1-b(i))^2/(cosd(gamma_prime(i))^2)) *x(i)^3;
    Cp_dx(i) = (pi^4/4) * s(i) * Cl*cosd(gamma_prime(i)) * ((1-b(i))^2/(cosd(gamma_prime(i))^2)) *x(i)^4;
    
    
end

C_T= trapz(x,Ct_dx);
C_P=trapz(x,Cp_dx);

T = C_T * rho * n^2 * d^4;

end

% FUNCTION TO GET CL AND CD IN THE BEMT part of the code from lab 3
function [Cl, Cd] = GETClCd(alpha_prime, Re, Re_values, AoA_vector, Cl_NACA, Cd_NACA)

[~, i] = min(abs(Re_values - Re));
Cl_vector = Cl_NACA(:,i);
Cd_vector = Cd_NACA(:,i);

if alpha_prime < AoA_vector(1)
    alpha_prime = AoA_vector(1);
elseif alpha_prime > AoA_vector(end)
    alpha_prime = AoA_vector(end);
end

Cl = interp1(AoA_vector, Cl_vector, alpha_prime, 'linear', 'extrap');
Cd = interp1(AoA_vector, Cd_vector, alpha_prime, 'linear', 'extrap');
end