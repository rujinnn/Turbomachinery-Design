%% TURBOPROP DESING
%% LAB 2

%% LAB GROUP J
%% MARIA SAN EMETERIO GONZALEZ  100471680
%% ANA BELEN ROMERA BAENA       100471676
%% WOOJIN LEE                   100576121
%% AIDAN JONES                  100578091

clear all
close all
clc

%% VARIABLES
% Chosen geometry 11x7 inches
geometryA = readmatrix('geometry11x7.txt');
r_R=geometryA(:,1);    % [m]
c_R=geometryA(:,2);    % [m]
beta=geometryA(:,3);   % [degrees]


% NACA0012 GIVEN BY PROFESSOR 
CD_NACA0012 = readmatrix('CDNACA0012.txt');
CL_NACA0012 = readmatrix('CLNACA0012.txt');

% DATA CHOSEN for 3 diferent RPM 3013, 4017 and 5018
data1 = readmatrix('data3013.txt');
data2 = readmatrix('data4017.txt');
data3 = readmatrix('data5018.txt');

% store the different values in their respective matrices creatin 17x3
% matricces

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


% rest of the variables 
D=0.2798;         % diammeter         [m]
R=D/2;            % radius            [m]
mu = 1.81*10^-5;  % dynamic viscosity [kg/(ms)]
rho = 1.225;      % air desnsity      [kg/m^3]
n=2;              % number of blades 

rpm=[3013,4017,5018]; %revolutions per minute

%% CALCULATIONS 

% k goes from 1 to 3, since we have 3 RPM
for k=1:3
    
    N=rpm(k)/60;      % rotational speed [RPS] 
    
    % gamma, alpha(AoA), V, and Re for each section of the blade depending
    % on J and r_R both values given in the data sets 
    
    % also calculates Cl, Cd,  dL, dD
    
    for i=1:length(J_data)
        for j=1:length(r_R)
            gamma(i,j) = atand(J_data(i,k)/(pi*r_R(j)));                  % [degrees]
            alpha(i,j) = beta(j)-gamma(i,j);                              % alpha along the blade [degrees]
            Vr(i,j)=sqrt(N^2*D^2*(J_data(i,k)^2+((2*pi*r_R(j)*R)/D)^2));  % resultant speed along the blade
            Re_Vr(i,j)= (Vr(i,j)*D*rho)/mu;                               % Re for the resultant speed
            
            [Cl, Cd] = ClCd(alpha(i,j), Re_Vr(i,j), CL_NACA0012, CD_NACA0012);  % function for the interpolation 
            
            
            dL(i,j) = Cl*(rho/2)*N^2*D^2*(J_data(i).^2+((2*pi*r_R(j)*R)/D).^2).*c_R(j).*R;
            dD(i,j) = Cd*(rho/2)*N^2*D^2*(J_data(i).^2+((2*pi*r_R(j)*R)/D).^2).*c_R(j).*R;
        end
    end
      
    
    dT = dL.*cosd(gamma)-dD.*sind(gamma);                          % dT = (dL cos(gamma) - dD sin(gamma))
    dP = (dL.*sind(gamma)+dD.*cosd(gamma)).*(D*N*pi*r_R)';       % dT = (dL cos(gamma) - dD sin(gamma))
    
    % intagrate using MATLAB trapezoidal function to obtain T and P
    
    % calculate Cp and Ct and use the results to compute the efficiency for
    % both the experimental and the BET results
        
    % THRUST
    T = n*trapz(r_R*R, dT, 2);
    CT_results = T./(rho*N^2*D^4);
    
    % POWER
    P = n*trapz(r_R*R, dP, 2);
    CP_results = P./(rho*N^3*D^5);
    
    % EFFICIENCY
    eta(:,k) = (CT_results./CP_results).*J_data(:,k);
    
 
    CT(k,:)=CT_results;
    CP(k,:)=CP_results;
end


%% PLOT BET RESULTS AGAINST EXPERIMENTAL ONES FOR COMPARISON 
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% CT
figure('Position', [100, 100, 1400, 450])
tit = {'$C_T~vs~J~for~n=3013~rpm$', '$C_T~vs~J~for~n=4017~rpm$', '$C_T~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), CT(i,:),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CT_data(:,i),'-r', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$C_T$')
    grid on
    title(tit{i})
    legend('BET', 'Experimental', 'Location', 'best')
    xlim([0 0.9])
    set(gca, 'FontSize', 12)
end
saveas(gcf, 'grafica_CT_vs_J.png')


%CP
figure('Position', [100, 100, 1400, 450])
tit = {'$C_P~vs~J~for~n=3013~rpm$', '$C_P~vs~J~for~n=4017~rpm$', '$C_P~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), CP(i,:),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), CP_data(:,i),'-r', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$C_P$')
    grid on
    title(tit{i})
    legend('BET', 'Experimental', 'Location', 'best')
    set(gca, 'FontSize', 12)
    
end
saveas(gcf, 'grafica_CP_vs_J.png')

%eta
figure('Position', [100, 100, 1400, 450])
tit = {'$eta~vs~J~for~n=3013~rpm$', '$eta~vs~J~for~n=4017~rpm$', '$eta~vs~J~for~n=5018~rpm$'};
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'compact')
for i = 1:3
    nexttile
    plot(J_data(:,i), eta(:,i),'-b', 'LineWidth', 1.5)
    hold on
    plot(J_data(:,i), eta_data(:,i),'-r', 'LineWidth', 1.5)
    xlabel('$J$')
    ylabel('$\eta$')
    grid on
    title(tit{i})
    legend('BET', 'Experimental', 'Location', 'best')
end
saveas(gcf, 'grafica_eta_vs_J.png')

%% PLOT BET RESULTR FOR DIFFERENT RPM

% CT
figure(4)
plot(J_data(:,1), CT(1,:),'-b', 'LineWidth', 1.5)
hold on
plot(J_data(:,2), CT(2,:),'-g', 'LineWidth', 1.5)
hold on
plot(J_data(:,3), CT(3,:),'-r', 'LineWidth', 1.5)
title('$C_T~vs~J~$')
legend('3013', '4017','5018','location', 'best')
grid on
xlabel('$J$')
ylabel('$C_T$')
saveas(figure(4),'Ct_vs_Jall.png');

%CP
figure (5)
plot(J_data(:,1), CP(1,:),'-b','LineWidth', 1.5)
hold on
plot(J_data(:,2), CP(2,:),'-g','LineWidth', 1.5)
hold on
plot(J_data(:,3), CP(3,:),'-r','LineWidth', 1.5)
xlabel('$J$')
ylabel('$C_P$')
grid on
title('$C_P~vs~J~$')
legend('3013', '4017','5018','location', 'best')
saveas(figure(5),'Cp_vs_Jall.png');

%eta
figure (6)
plot(J_data(:,1), eta(:,1),'-b', 'LineWidth', 1.5)
hold on
plot(J_data(:,2), eta(:,2),'-g', 'LineWidth', 1.5)
hold on
plot(J_data(:,3), eta(:,3),'-r', 'LineWidth', 1.5)
xlabel('$J$')
ylabel('$eta$')
grid on
title('$\eta~vs~J~$')
legend('3013', '4017','5018', 'location', 'best')
saveas(figure(6),'eta_vs_Jall.png');

%% FUNCTION 
function [Cl, Cd] = ClCd(alpha, Re,CL_NACA0012, CD_NACA0012)

% obtaining the angle of attack from the NACA0012 files (frist column)
alpha_CDNACA = CD_NACA0012(:,1);
alpha_CLNACA = CL_NACA0012(:,1);

% Re number from the NACA0012  data file 
Re_NA0012 = [160000, 360000, 700000, 100000, 2000000, 5000000];

% using 2 iterpolations 
% first interpolating Re keeping alpha constant
% and then interpolating alpha keeping Re constant

% interpolating Re with alpha constant
for i=1:length(alpha_CDNACA)
    Cd_AoA_cte(i,:) = interp1(Re_NA0012, CD_NACA0012(i,2:end), Re, 'makima');
    Cl_AoA_cte(i,:) = interp1(Re_NA0012, CL_NACA0012(i,2:end), Re, 'makima');
end

% Interpolating alpha with Re constant
Cd = interp1(alpha_CDNACA', Cd_AoA_cte, alpha, 'makima');
Cl = interp1(alpha_CLNACA', Cl_AoA_cte, alpha, 'makima');

end
