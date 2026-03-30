%%Inputs 
P = 101.33; % kPa

% Liquid composition sweep
x1_vals = linspace(0, 1, 50); 

% NRTL parameters
b12 = 631.05; b21 = 1197.41; alpha = 0.5343;
parameters = [b12,b21,alpha];

% Antoine parameters (grouped in two vectors)
ant1 = [14.3145, 2756.22, 228.06];   % component 1
ant2 = [16.3872, 3885.70, 230.170];  % component 2

%%%Equilibrium plot using NRTL
y1star = equilibrium(x1_vals, P, ant1, ant2, parameters);

%%%Mccabe Inputs
R = 3;
xd = .9;
xb = .1;
xz = .5;
q = .5;

% Liquid mole fraction Experimental Data
x = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1];

% Vapor mole fraction Experimental Data
y = [0, 0.6381, 0.7301, 0.7716, 0.7916, 0.8124, 0.8269, 0.8387, 0.8532, 0.8712, 0.895, 0.9335, 0.9627, 1];

plot(x1_vals,equilibrium(x1_vals, P, ant1, ant2, parameters))
hold on
plot(linspace(0,1,10),linspace(0,1,10))
plot(x,y)

%Pretty Good Fit! MSE is about 1%

%%Distillation Logic%%
% Rectifying Equation
y_rect = @(x) (R/(R+1))*x + xd/(R+1);

% q-line
if q == 1
    x_q = xz; % vertical line
    plot([x_q x_q], [0 1], 'm--', 'LineWidth',2);
    y_intersect = y_rect(x_q);
else
    y_q = @(x) (q/(q-1))*x - xz/(q-1);
    plot(x1_vals, y_q(x1_vals), 'm--', 'LineWidth',2);
    
    % Intersection
    x_q = fsolve(@(x) y_rect(x) - y_q(x), xz);
    y_intersect = y_rect(x_q);
end

% Stripping Equation
m_strip = (y_intersect - xb)/(x_q - xb);
y_strip = @(x) m_strip*(x - xb) + xb;

%Plotting the operating lines
plot([xd,x_q],[xd,y_intersect])
plot([xb,x_q], [xb, y_intersect], 'g', 'LineWidth',2);

%Stages
stage = 0;
x_stage = xd;

y_eq = @(x) interp1(x1_vals, y1star, x, 'pchip');

while x_stage > xb
    % Step 1: Horizontal to equilibrium curve
    if x_stage >= x_q
        y_stage = y_rect(x_stage); % rectifying
    else
        y_stage = y_strip(x_stage); % stripping
    end
    f = @(x) interp1(x1_vals, y1star, x) - y_stage;
    x_stage_star = fzero(f, [0,x_stage]); % initial guess
    
    % Plot horizontal
    plot([x_stage, x_stage_star], [y_stage, y_stage], 'm--');
    
    % Step 2: Vertical to operating line
    if x_stage_star >= x_q
        y_next = y_rect(x_stage_star); % rectifying
        else
            y_next = y_strip(x_stage_star); % stripping
    end
    
    %Plot Vertical 
    plot([x_stage_star, x_stage_star],[y_stage, y_next], 'm--')
    x_stage_previous = x_stage;
    x_stage = x_stage_star;
    stage = stage + 1;
    % Distillation probably not appropriate at this point >:(
    if stage > 50
        break
    end
end

stage = stage - (xb-x_stage)/(x_stage_previous - x_stage); %Fractional Stage Calculation

xlabel("X1")
ylabel("y1")
title("McCabe Thiele Analysis")
legend("NRTL Model", "1 to 1 Line", "Experimental Data", "Q line", "rectifying Line", "Stripping Line", 'Location','southeast')
hold off

%Sensitivity Analysis Between R and N

%Minimum Rectifying Line
%Should work iwth both a dataset and just plain old NRTL paramters
%y1star is the corresponding equilibrium values calculated using NRTL if
%done that way
R_min = NRTL_minimum_rect(x1_vals, y1star, xd, xz,q);
%Minimum Stripping Line
B_min = NRTL_minimum_strip(x1_vals, y1star, xb, xz,q);

R1 = 1.1*R_min;
R2 = 3*R_min;
B1 = 1.1*B_min;
B2 = 3*B_min;

points = 50;
R = linspace(R1, R2, points);
B = linspace(B1, B2, points);

stages = linspace(1,10,points);

%Specifying Reflux Ratio
for i = 1:length(R)
    y_rect = @(x) (R(i)/(R(i)+1))*x + xd/(R(i)+1);
    % q-line
    if q == 1
        x_q = xz; % vertical line
        y_intersect = y_rect(x_q);
    else
        y_q = @(x) (q/(q-1))*x - xz/(q-1);
        % Intersection
        x_q = fsolve(@(x) y_rect(x) - y_q(x), xz);
        y_intersect = y_rect(x_q);
    end

    % Stripping Equation
    m_strip = (y_intersect - xb)/(x_q - xb);
    y_strip = @(x) m_strip*(x - xb) + xb;
    
    %Defining initial Conditions and Equilibrium Cruve for Stage counting
    stages(i) = 0;
    x_stage = xd;

    y_eq = @(x) interp1(x1_vals, y1star, x, 'pchip');
    while x_stage > xb
        % Step 1: Horizontal to equilibrium curve
        if x_stage >= x_q
            y_stage = y_rect(x_stage); % rectifying
        else
            y_stage = y_strip(x_stage); % stripping
        end
        f = @(x) interp1(x1_vals, y1star, x) - y_stage;
        x_stage_star = fzero(f, [0,x_stage]); % initial guess
        
        % Step 2: Vertical to operating line
        if x_stage_star >= x_q
            y_next = y_rect(x_stage_star); % rectifying
            else
                y_next = y_strip(x_stage_star); % stripping
        end
        x_stage_previous = x_stage;
        x_stage = x_stage_star;
        stages(i) = stages(i) + 1;
        % Distillation probably not appropriate at this point >:(
        if stage > 50
            break
        end
    end
    stages(i) = stages(i) - (xb-x_stage)/(x_stage_previous - x_stage);
end
hold off

%Plotting the senstitivity analysis
%{
plot(R,stages)
hold on
xlabel("Reflux Ratio")
ylabel("Number of Stages")
title("Sensitivity Analysis")
hold off
%}

%Function that I plug into the solver
function F = nrtl_bubble(T, x1, x2, P, b12, b21, alpha, ant1, ant2, R)
    % NRTL activity coefficients Calculations
    tau12 = b12/(R*T);
    tau21 = b21/(R*T);

    G12 = exp(-alpha*tau12);
    G21 = exp(-alpha*tau21);

    gamma1 = exp(x2^2*((tau21*(G21/(x1+x2*G21))^2) + (tau12*G12/(x2+x1*G12)^2)));
    gamma2 = exp(x1^2*((tau12*(G12/(x2+x1*G12))^2) + (tau21*G21/(x1+x2*G21))^2));

    P1sat = exp(ant1(1) - ant1(2)/(T - 273.15 + ant1(3)));
    P2sat = exp(ant2(1) - ant2(2)/(T -273.15 + ant2(3)));

    % Bubble point Equation
    F = P - (x1*gamma1*P1sat + x2*gamma2*P2sat);
end
function equilibrium = equilibrium(x1_vals, P, ant1, ant2, parameters)
    %%Initialization
    x2_vals = 1-x1_vals;
    equilibrium = x1_vals;
    T_vals = x1_vals;
    
    % Gas constant in cal/mol-K
    R = 8.314/4.184;
    
    %Pure Component Boiling Points
    T_b1 = fsolve(@(T) P - exp(ant1(1) - ant1(2)/(T -273.15+ ant1(3))), 350);
    T_b2 = fsolve(@(T) P - exp(ant2(1) - ant2(2)/(T -273.15+ ant2(3))), 350);
    
    %paramters
    b12 = parameters(1);
    b21 = parameters(2);
    alpha = parameters(3);

    %%%%%For loop Logic for Equilibrium
    for i = 1:length(x1_vals)
        x1 = x1_vals(i);
        x2 = x2_vals(i);
        
        % Initial guess: interpolate between pure component boiling points
        T_guess = T_b1 + x1*(T_b2 - T_b1);
    
        % Solve bubble point temperature using a solver
        T_vals(i) = fsolve(@(T) nrtl_bubble(T, x1, x2, P, b12, b21, alpha, ant1, ant2, R), T_guess);
    
        % Solve for saturation pressure at given bubble point
        P1sat = exp(ant1(1) - ant1(2)/(T_vals(i) - 273.15 + ant1(3)));
        P2sat = exp(ant2(1) - ant2(2)/(T_vals(i) - 273.15 + ant2(3)));
        
        %Solve for paramter tau
        tau12 = b12/(R*T_vals(i));
        tau21 = b21/(R*T_vals(i));
        
        %Solve for paramter G
        G12 = exp(-alpha*tau12);
        G21 = exp(-alpha*tau21);
        
        %Solve for activity coefficient
        gamma1 = exp(x2^2*((tau21*(G21/(x1+x2*G21))^2) + (tau12*G12/(x2+x1*G12)^2)));
        gamma2 = exp(x1^2*((tau12*(G12/(x2+x1*G12))^2) + (tau21*G21/(x1+x2*G21))^2));
        
        %Relatioship between partial and total Pressure
        equilibrium(i) = x1*gamma1*P1sat / P;
    end
    %%%%%%%%
end

function min_reflux = NRTL_minimum_rect(x1_vals, y1star, xd, xz,q)
    %interpolated function
    y_eq = @(x) interp1(x1_vals, y1star, x, 'pchip');

    if q == 1
        x_q = xz; % vertical line
        y_coord = y_eq(x_q);
        min_reflux = (xd - y_coord) / (y_coord - x_q);
    else
        y_q = @(x) (q/(q-1))*x - xz/(q-1);
        x1star = fsolve(@(x) y_eq(x) - y_q(x), xz);
        y_coord = y_eq(x1star);
        min_reflux = (xd - y_coord) / (y_coord - x1star);
    end
    
end

function min_strip = NRTL_minimum_strip(x1_vals, y1star, xb, xz,q)
    %interpolated function
    y_eq = @(x) interp1(x1_vals, y1star, x, 'pchip');

    if q == 1
        x_q = xz; % vertical line
        y_coord = y_eq(x_q);
        min_strip = (y_coord - xb) / (x_q - xb);
    else
        y_q = @(x) (q/(q-1))*x - xz/(q-1);
        x1star = fsolve(@(x) y_eq(x) - y_q(x), xz);
        y_coord = y_eq(x1star);
        min_strip= (y_coord - xb) / (x1star - xb);
    end
end

function min_stages = NRTL_minimum_stages(x1_vals, y1star, xd, xb)
    x_stage = xd;
    
    min_stages= 0;
    while x_stage > xb
    % Step 1: Horizontal to equilibrium curve
    y_stage = x_stage;
    f = @(x) interp1(x1_vals, y1star, x) - y_stage;
    x_stage_star = fzero(f, [0,x_stage]); % initial guess
    
    
    % Plot horizontal
    plot([x_stage, x_stage_star], [y_stage, y_stage]);
    
    % Step 2: Vertical to operating line
    
    %Plot Vertical 
    plot([x_stage_star, x_stage_star],[y_stage, x_stage_star])

    x_stage = x_stage_star;
    min_stages = min_stages + 1;
    % Distillation probably not appropriate at this point >:(
    if min_stages > 50
        break
    end
    end
end
function MSE = MSE(x, y, P, ant1, ant2, parameters)
    model_y = x;
    differences = 0;
    for i = 1:length(x)
        model_y(i) = equilibrium(x(i), P, ant1, ant2, parameters);
        differences = differences + (y(i)*100 - model_y(i)*100)^2;
    end
    MSE = differences/length(x);
end

