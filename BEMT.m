%% TURBOPROP DESING
%% LAB 3: PROPELLER PEROFORMANCES WITH BEMT

%% LAB GROUP J
%% MARIA SAN EMETERIO GONZALEZ  100471680
%% ANA BELEN ROMERA BAENA       100471676
%% WOOJIN LEE                   100576121
%% AIDAN JONES

clear all
close all
clc

%% GRAPHS in latex format for the report
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% VARIABLES
%% GEOMETRY
% Chosen geometry 11x7 inches
% same as in lab 2
geometryA = readmatrix('geometry11x7.txt');
r_R=geometryA(:,1);    % [m]
c_R=geometryA(:,2);    % [m]
beta=geometryA(:,3);   % [degrees]

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

%% RPM 
% DATA CHOSEN for 3 diferent RPM 3013, 4017 and 5018
data1 = readmatrix('data3013.txt');
data2 = readmatrix('data4017.txt');
data3 = readmatrix('data5018.txt');

% extract data from the previous matrices
J_data(:,1) = data1(:,1);  
J_data(:,2) = data2(:,1);
J_data(:,3) = data3(:,1);
CT_data(:,1) = data1(:, 2); 
CT_data(:,2) = data2(:, 2);
CT_data(:,3) = data3(:, 2);
CP_data(:,1) = data1(:, 3); 
CP_data(:,2) = data2(:, 3);
CP_data(:,3) = data3(:, 3);
eta_data(:,1) = data1(:, 4); 
eta_data(:,2) = data2(:, 4);
eta_data(:,3) = data3(:, 4);
data_vector = {data1, data2, data3};

%% rest of the variables 
% From the propeller
D=0.2798;         % diammeter         [m]
R=D/2;            % radius            [m]
mu = 1.81e-5;  % dynamic viscosity [kg/(ms)]
rho = 1.225;      % air desnsity      [kg/m^3]
n=2;              % number of blades 

RPM_vector = [3013, 4017, 5018]; 
num_J = length(data1(:,1));   % 17
num_r = length(r_R);          % 18
num_RPM = length(RPM_vector); % 3

max_iter = 1e8; % max nummber of iterations

% Inicialize matrices for the results
a_f = zeros(num_r, num_J, num_RPM);
b_f = zeros(num_r, num_J, num_RPM);
gamma_prime_f = zeros(num_r, num_J, num_RPM);
alpha_prime_f = zeros(num_r, num_J, num_RPM);
Re_f = zeros(num_r, num_J, num_RPM);
Cl_f = zeros(num_r, num_J, num_RPM);
Cd_f = zeros(num_r, num_J, num_RPM);
Ve_f = zeros(num_r, num_J, num_RPM);

%% CALCULATIONS
% Solidity and x 
x = r_R;
s = (n*c_R)./(2*pi*r_R);

for k = 1:num_RPM % from 1 to 3 because we're testing 3 velocities 
    rpm= RPM_vector(k); % revolutions per minute
    N=rpm/60;           % revolutions per second
    w=N*2*pi;           % omega in          [rad/s]
    
    J = data_vector{k}(:,1); 
    
    for i = 1:length(J)
        
        % initial values for gamma a and b 
        gamma_prime = 10;
        a = 0.001;
        b = 0.0001;
        % tolerance 
        tol = 1e-8;
        
        for j = 1:length(x)
            iter = 0;
            converged = false; % inizialice converged to false so it enter in the while loop
            
            while ~converged && iter <= max_iter
                iter = iter + 1; 
                
                % STEP 1
                alpha_prime = beta(j) - gamma_prime;
                
                Ve = sqrt( ( J(i) * N * D )^2 * (1+a)^2 + ( w * x(j) * R )^2 * (1-b)^2 ); % effective velocity
                Re = (rho * Ve * c_R(j) * R) / mu; % Re number 
                
                % STEP 2
                [Cl, Cd] = GETClCd(alpha_prime, Re, Re_NACA0012, alpha_NACA, Cl_NACA, Cd_NACA);  % calling the function for calculating cl and d
                
                % STEP 3
                % calculate a and b for next iteration
                a_i = (s(j) / (4 * sind(gamma_prime)^2) * (Cl * cosd(gamma_prime) - Cd * sind(gamma_prime)));
                % a_i=a_cal/1+a_cal
                b_i = (s(j)/ (4 * sind(gamma_prime) * cosd(gamma_prime)) * (Cl * sind(gamma_prime) + Cd * cosd(gamma_prime)));
                % b_i=b_cal/1-b_cal
                
                a_cal = a_i/(1-a_i);
                b_cal = b_i/(b_i+1);
                
                % STEP 4
                J_cal = (pi*x(j))*((1-b)/(1+a))*tand(gamma_prime); 
                
                % STEP 5
                if abs(J_cal-J(i))<tol
                    converged = true;
                else
                    % Change the values for a, b and gamma so it enter the
                    % loop again until J_obj is reach
                    a = a_cal; 
                    b = b_cal;
                    gamma_prime = gamma_prime + 0.1 * (J(i)-J_cal);
                end
            end
            
            if converged          
               % store the calculated values in matrices
                a_f(j,i,k) = a;
                b_f(j,i,k) = b;
                gamma_prime_f(j,i,k) = gamma_prime;
                alpha_prime_f(j,i,k) = alpha_prime;
                Re_f(j,i,k) = Re;
                Cl_f(j,i,k) = Cl;
                Cd_f(j,i,k) = Cd;
                Ve_f(j,i,k) = Ve;
            else
                % NO CONVERGED 
                fprintf('did not converged'); % display message that the iteration did not converged
            end
            
            % STEP 6
            % dCT/dx and dCP/dx
            dC_T(j,i,k) = (pi^3/4) * s(j) * (Cl_f(j,i,k) * cosd(gamma_prime_f(j,i,k)) - Cd_f(j,i,k) * sind(gamma_prime_f(j,i,k)))*...
                ((1-b_f(j,i,k))^2/cosd(gamma_prime_f(j,i,k))^2) * x(j)^3;
           
            dC_P (j,i,k)= ((pi^4)/4) * s(j) * (Cl_f(j,i,k) * sind(gamma_prime_f(j,i,k)) + Cd_f(j,i,k) * cosd(gamma_prime_f(j,i,k)))*...
                ((1-b_f(j,i,k))^2/cosd(gamma_prime_f(j,i,k))^2) * x(j)^4;
        end
                
        % STEP 8
        % Integrate over x using the trapezoidal function from MATLAB  
        C_T(i,k)= trapz(x,dC_T(:,i,k));
        C_P(i,k)=trapz(x,dC_P(:,i,k));
        eta(i,k)=(C_T(i,k)/C_P(i,k)) * J_data(i,k);
        
    end
end

%% VALUES FROM LAB 2
% extracted directly from last laboratory
CT_BET=readmatrix('CT_BET.txt');
CP_BET=readmatrix('CP_BET.txt');
eta_BET=readmatrix('eta_BET.txt');

%% PLOTS 
% plot BEMT  all velocities 
figure()
plot(J_data(:,1),C_T(:,1),'-b', 'LineWidth', 1.5)
hold on
plot(J_data(:,2),C_T(:,2),'-r', 'LineWidth', 1.5)
hold on
plot(J_data(:,3),C_T(:,3),'-g', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$C_T$')
grid on
title('$C_T~vs~J$')
legend('3013', '4017','5018', 'Location', 'best')
set(gca, 'FontSize', 12)
saveas(gcf, 'grafica_CT_vs_J.png')

figure()
plot(J_data(:,1),C_P(:,1),'-b', 'LineWidth', 1.5)
hold on
plot(J_data(:,2),C_P(:,2),'-r', 'LineWidth', 1.5)
hold on
plot(J_data(:,3),C_P(:,3),'-g', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$C_P$')
grid on
title('$C_P~vs~J$')
legend('3013', '4017','5018', 'Location', 'best')
set(gca, 'FontSize', 12)
saveas(gcf, 'grafica_CP_vs_J.png')

figure()
plot(J_data(:,1),eta(:,1),'-b', 'LineWidth', 1.5)
hold on
plot(J_data(:,2),eta(:,2),'-r', 'LineWidth', 1.5)
hold on
plot(J_data(:,3),eta(:,3),'-g', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$\eta$')
grid on
title('$\eta~vs~J$')
legend('3013', '4017','5018', 'Location', 'best')
set(gca, 'FontSize', 12)
saveas(gcf, 'grafica_eta_vs_J.png')

% plot BEMT  with BET and the experimental data for all RPM 
figure('Position', [100, 100, 1400, 450])
tit = {'$C_T~vs~J~for~n=3013~rpm$', '$C_T~vs~J~for~n=4017~rpm$', '$C_T~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), C_T(:,i),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CT_data(:,i),'-r', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CT_BET(i,:),'-g', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$C_T$')
    grid on
    title(tit{i})
    legend('BEMT', 'Experimental','BET', 'Location', 'best')
    set(gca, 'FontSize', 12)
end
saveas(gcf, 'grafica_CT_vs_J_all.png')

figure('Position', [100, 100, 1400, 450])
tit = {'$C_P~vs~J~for~n=3013~rpm$', '$C_P~vs~J~for~n=4017~rpm$', '$C_P~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), C_P(:,i),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CP_data(:,i),'-r', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CP_BET(i,:),'-g', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$C_P$')
    grid on
    title(tit{i})
    legend('BEMT', 'Experimental','BET', 'Location', 'best')
    set(gca, 'FontSize', 12)
end
saveas(gcf, 'grafica_CP_vs_J_all.png')

figure('Position', [100, 100, 1400, 450])
tit = {'$\eta~vs~J~for~n=3013~rpm$', '$\eta~vs~J~for~n=4017~rpm$', '$\eta~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), eta(:,i),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), eta_data(:,i),'-r', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), eta_BET(:,i),'-g', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$\eta$')
    grid on
    title(tit{i})
    legend('BEMT', 'Experimental','BET', 'Location', 'best')
    set(gca, 'FontSize', 12)
end
saveas(gcf, 'grafica_eta_vs_J_all.png')






%% FUNCTION
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
