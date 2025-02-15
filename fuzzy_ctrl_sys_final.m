% Bryan Esparza, Gian Fajardo, Ejmin Gharibian
% Dr. Jia
% ECE 581
% 
% Final Project - Advanced Cruise Control

clear, clc, close all
format compact, format short



% OUTLINE: ---------------------------------------------------------------
% figure management, then
% 
% the fuzzy systems of:
%
%   the braking distance FIS, featuring:
%     Velocity_in input membership functions,
%     Braking_Distance output membership functions, 
%
%     both of their MF graphs, and
%     the input/output relationship.
%
% 
%   the fuzzy cruise controller, featuring:
%     Position_error input membership functions,
%     Velocity_error input membership functions,
%     velocity_output output membership functions,
%     accel_out output membership functions,
% 
%     the input and output MF graphs, and
%     the surface graphs of:
%       the velocity_output output, and
%       the accel_out output. Then,
%       
% 
% the state of the front car, then
% the state space model layout of the model ACC_car, then
% the additional animation layer, then
% the main control loop, where:
%   all of the states of ACC_car are processed,
%   the animation is displayed onto a figure, and
%   the states of everything are saved onto graphs using graph_array.
%
% Finally, all the state-related graphs (a, v, x) are displayed




% figure management ------------------------------------------------------
% this helps us manage the many figures of this MATLAB script

% If true, each sets of figures will load sequentially. If false, then
%   it'll skip to the results.
tour_each_subsys  = true;
window_split      = true;

% Uncomment which screen dimensions are right for you!
screen_dims     = [2880 1864];  % 2024 MacBook Air 15' screen resolution
% screen_dims     = [1920 1080];  % typical Windows laptop screen dimensions




num_figures_x_y = [5 3];        % num of figures along x and y dir of screen

if window_split == true
  screen_positions_x = linspace( 0, screen_dims(1), num_figures_x_y(1)+1 );
  screen_positions_y = linspace( 0, screen_dims(2), num_figures_x_y(2)+1 );
else
  screen_positions_x = (screen_dims(1)/2) .* ones(num_figures_x_y(1));
  screen_positions_y = (screen_dims(2)/2) .* ones(num_figures_x_y(2));
end

% this runs alongside each figure
sp = @(x,y)  [screen_positions_x(x)  screen_positions_y(y)];











% fuzzy systems ----------------------------------------------------------

% Important for setting the limits of the universe of discourse after
%   each membership function triangle widths are made.
min_max = @(vector) [min(vector)  max(vector)];

% goes from 1 to 2; no rain: 1, max rain: 2
rain_intensity  = 1;      % multiplier

bd_update_rate  = 10;     % number of samples
is_it_raining   = false;  % term which activates raining simulation

% restricts rain_intensity 'x' between a and b
limiter = @(x, a,b)  min( max(a, x), b );


% braking distance FIS
bd_fis  = mamfis(                         ...
    'NumInputs',  1, 'NumInputMFs',   4,  ...
    'NumOutputs', 1, 'NumOutputMFs',  4,  ...
    'AddRule', 'none'                     ...
);

% for your reference, the following spreadsheet is the table used to link 
%   the two sets of numbers together. 
% https://docs.google.com/spreadsheets/d/1iIUxjGq-N72JQ1pvdXxdex4cZxI1lnz1IOF_hztId_Y/edit?usp=sharing


% Velocity_in membership functions ------------------------------------
VI_peaks        = [0.0  13.4  22.4  31.3];
distances       = VI_peaks(2:end) - VI_peaks(1:(end-1));

mf1_parameters  = [(-distances(1) + VI_peaks(1))  VI_peaks(1:2)];
mf2_parameters  = VI_peaks(1:3);
mf3_parameters  = VI_peaks(2:4);
mf4_parameters  = [VI_peaks(3:4)  (VI_peaks(4) + distances(end))];


VI_universe     = min_max( [mf1_parameters  mf4_parameters] )

bd_fis.Inputs(1).Name = 'Velocity_in';
bd_fis.Inputs(1).Range = VI_universe;

bd_fis.Inputs(1).MembershipFunctions(1).Name       = 'Zero';
bd_fis.Inputs(1).MembershipFunctions(1).Type       = 'trimf'; 
bd_fis.Inputs(1).MembershipFunctions(1).Parameters = mf1_parameters;

bd_fis.Inputs(1).MembershipFunctions(2).Name       = 'Low';
bd_fis.Inputs(1).MembershipFunctions(2).Type       = 'trimf';
bd_fis.Inputs(1).MembershipFunctions(2).Parameters = mf2_parameters;

bd_fis.Inputs(1).MembershipFunctions(3).Name       = 'Medium';
bd_fis.Inputs(1).MembershipFunctions(3).Type       = 'trimf';
bd_fis.Inputs(1).MembershipFunctions(3).Parameters = mf3_parameters;

bd_fis.Inputs(1).MembershipFunctions(4).Name       = 'High';
bd_fis.Inputs(1).MembershipFunctions(4).Type       = 'trimf';
bd_fis.Inputs(1).MembershipFunctions(4).Parameters = mf4_parameters;



% Braking_Distance membership functions ------------------------------------
BD_peaks        = [1.5  27.1  60.4  106.1];
distances       = BD_peaks(2:end) - BD_peaks(1:(end-1));

mf1_parameters  = [(-distances(1) + BD_peaks(1))  BD_peaks(1:2)];
mf2_parameters  = BD_peaks(1:3);
mf3_parameters  = BD_peaks(2:4);
mf4_parameters  = [BD_peaks(3:4)  (BD_peaks(4) + distances(end))];


% my anonymous function. see above
BD_universe     = min_max( [mf1_parameters  mf4_parameters] )

bd_fis.Outputs(1).Name  = 'Braking_Distance';
bd_fis.Outputs(1).Range = BD_universe;

bd_fis.Outputs(1).MembershipFunctions(1).Name       = 'Close';
bd_fis.Outputs(1).MembershipFunctions(1).Type       = 'trimf'; 
bd_fis.Outputs(1).MembershipFunctions(1).Parameters = mf1_parameters;

bd_fis.Outputs(1).MembershipFunctions(2).Name       = 'Slight';
bd_fis.Outputs(1).MembershipFunctions(2).Type       = 'trimf';
bd_fis.Outputs(1).MembershipFunctions(2).Parameters = mf2_parameters;

bd_fis.Outputs(1).MembershipFunctions(3).Name       = 'Medium';
bd_fis.Outputs(1).MembershipFunctions(3).Type       = 'trimf';
bd_fis.Outputs(1).MembershipFunctions(3).Parameters = mf3_parameters;

bd_fis.Outputs(1).MembershipFunctions(4).Name       = 'Far';
bd_fis.Outputs(1).MembershipFunctions(4).Type       = 'trimf';
bd_fis.Outputs(1).MembershipFunctions(4).Parameters = mf4_parameters;

ruleset_braking_distance = [
  "if Velocity_in is Zero then Braking_Distance is Close"
  "if Velocity_in is Low then Braking_Distance is Slight"
  "if Velocity_in is Medium then Braking_Distance is Medium"
  "if Velocity_in is High then Braking_Distance is Far"
];

bd_fis = addRule(bd_fis, ruleset_braking_distance);

clear ruleset_braking_distance


f1 = figure; f1.Position(1:2) = sp(1,1);
  plotmf(bd_fis, 'input', 1, 250);
  title 'Velocity_in Membership Functions'

f1 = figure; f1.Position(1:2) = sp(1,2);
  plotmf(bd_fis, 'output', 1, 250);
  title 'Braking Distance Membership Functions'


% testing the output to see if the in-between values are being interpolated
velocity_inputs = ( VI_peaks(1):VI_peaks(end) )';
BD_outputs      = evalfis(bd_fis, velocity_inputs);


f1 = figure; f1.Position(1:2) = sp(2,1);
  scatter(VI_peaks, BD_peaks, 'filled'), hold on, grid on
  plot(velocity_inputs, BD_outputs, LineWidth=2)
    title 'Braking Distance as a function of Velocity'
    xlabel 'Velocity [m/s]', ylabel 'Braking Distance [m]'

clear velocity_inputs BD_outputs

pause_tour(tour_each_subsys, "braking distance FIS is done.")





acc_fis  = mamfis(                        ...
    'NumInputs',  2, 'NumInputMFs',   3,  ...
    'NumOutputs', 2, 'NumOutputMFs',  3,  ...
    'AddRule', 'none'                     ...
);

% Position_Error membership functions ------------------------------------
PE_width        = 500;   % m; determines the acceptable input range

mf1_parameters  = [-1/2  0  1/2].*PE_width;
mf2_parameters  = [-1/2  0  1/2].*PE_width + PE_width/2;
mf3_parameters  = [-1/2  0  1/2].*PE_width + PE_width;

% my anonymous function. see above
PE_universe     = min_max( [mf1_parameters  mf3_parameters] )


acc_fis.Inputs(1).Name = 'Position_Error';
acc_fis.Inputs(1).Range = PE_universe;

acc_fis.Inputs(1).MembershipFunctions(1).Name       = 'Zero';
acc_fis.Inputs(1).MembershipFunctions(1).Type       = 'trimf'; 
acc_fis.Inputs(1).MembershipFunctions(1).Parameters = mf1_parameters;

acc_fis.Inputs(1).MembershipFunctions(2).Name       = 'Moderate';
acc_fis.Inputs(1).MembershipFunctions(2).Type       = 'trimf';
acc_fis.Inputs(1).MembershipFunctions(2).Parameters = mf2_parameters;

acc_fis.Inputs(1).MembershipFunctions(3).Name       = 'High';
acc_fis.Inputs(1).MembershipFunctions(3).Type       = 'trimf';
acc_fis.Inputs(1).MembershipFunctions(3).Parameters = mf3_parameters;


f1 = figure; f1.Position(1:2) = sp(1,1);
  plotmf(acc_fis, 'input', 1, 250);
  title 'Position Error Membership Functions'



% Vel_Error membership functions -----------------------------------
VE_width        = 100;   % m/s; determines the acceptable input range

mf1_parameters  = [-1/2  0  1/2].*VE_width;
mf2_parameters  = [-1/2  0  1/2].*VE_width + VE_width/2;
mf3_parameters  = [-1/2  0  1/2].*VE_width + VE_width;


VE_universe     = min_max( [mf1_parameters  mf3_parameters] )


acc_fis.Inputs(2).Name = 'Vel_Error';
acc_fis.Inputs(2).Range = VE_universe;

acc_fis.Inputs(2).MembershipFunctions(1).Name       = 'Zero';
acc_fis.Inputs(2).MembershipFunctions(1).Type       = 'trimf'; 
acc_fis.Inputs(2).MembershipFunctions(1).Parameters = mf1_parameters;

acc_fis.Inputs(2).MembershipFunctions(2).Name       = 'Moderate';
acc_fis.Inputs(2).MembershipFunctions(2).Type       = 'trimf';
acc_fis.Inputs(2).MembershipFunctions(2).Parameters = mf2_parameters;

acc_fis.Inputs(2).MembershipFunctions(3).Name       = 'High';
acc_fis.Inputs(2).MembershipFunctions(3).Type       = 'trimf';
acc_fis.Inputs(2).MembershipFunctions(3).Parameters = mf3_parameters;


f1 = figure; f1.Position(1:2) = sp(1,2);
  plotmf(acc_fis, 'input', 2, 250);
  title 'Velocity Error Membership Functions'



% Output velocity membership functions -----------------------------------
VO_width        = 50;

mf1_parameters  = [-1/2  0  1/2].*VO_width;
mf2_parameters  = [-1/2  0  1/2].*VO_width + VO_width/2;
mf3_parameters  = [-1/2  0  1/2].*VO_width + VO_width;

% my anonymous function. see above
VO_universe     = min_max( [mf1_parameters  mf3_parameters] )


acc_fis.Outputs(1).Name = 'Change_Velocity';
acc_fis.Outputs(1).Range = VO_universe;

acc_fis.Outputs(1).MembershipFunctions(1).Name       = 'Zero';
acc_fis.Outputs(1).MembershipFunctions(1).Type       = 'trimf'; 
acc_fis.Outputs(1).MembershipFunctions(1).Parameters = mf1_parameters;

acc_fis.Outputs(1).MembershipFunctions(2).Name       = 'Moderate';
acc_fis.Outputs(1).MembershipFunctions(2).Type       = 'trimf';
acc_fis.Outputs(1).MembershipFunctions(2).Parameters = mf2_parameters;

acc_fis.Outputs(1).MembershipFunctions(3).Name       = 'High';
acc_fis.Outputs(1).MembershipFunctions(3).Type       = 'trimf';
acc_fis.Outputs(1).MembershipFunctions(3).Parameters = mf3_parameters;

f1 = figure; f1.Position(1:2) = sp(2,1);
  plotmf(acc_fis, 'output', 1, 250);
  title 'Change-Velocity Membership Functions'


% Output Acceleration membership functions -------------------------------
AO_width        = 4;

mf1_parameters  = [-1/2  0  1/2].*AO_width;
mf2_parameters  = [-1/2  0  1/2].*AO_width + AO_width/2;
mf3_parameters  = [-1/2  0  1/2].*AO_width + AO_width;

% my anonymous function. see above
AO_universe     = min_max( [mf1_parameters  mf3_parameters] )


acc_fis.Outputs(2).Name = 'Accel_Out';
acc_fis.Outputs(2).Range = AO_universe;

acc_fis.Outputs(2).MembershipFunctions(1).Name       = 'Zero';
acc_fis.Outputs(2).MembershipFunctions(1).Type       = 'trimf';
acc_fis.Outputs(2).MembershipFunctions(1).Parameters = mf1_parameters;

acc_fis.Outputs(2).MembershipFunctions(2).Name       = 'Moderate';
acc_fis.Outputs(2).MembershipFunctions(2).Type       = 'trimf';
acc_fis.Outputs(2).MembershipFunctions(2).Parameters = mf2_parameters;

acc_fis.Outputs(2).MembershipFunctions(3).Name       = 'High';
acc_fis.Outputs(2).MembershipFunctions(3).Type       = 'trimf';
acc_fis.Outputs(2).MembershipFunctions(3).Parameters = mf3_parameters;

f1 = figure; f1.Position(1:2) = sp(2,2);
  plotmf(acc_fis, 'output', 2, 250);
  title 'Accel-Out Membership Functions'


% Ruleset additions ------------------------------

ruleset_velocity = [
  "if Position_Error is Zero and Vel_Error is Zero then Change_Velocity is Zero"
  "if Position_Error is Zero and Vel_Error is not Zero then Change_Velocity is Moderate"
  "if Position_Error is not Zero and Vel_Error is Zero then Change_Velocity is Moderate"
  "if Position_Error is Moderate and Vel_Error is Moderate then Change_Velocity is Moderate"
  "if Position_Error is High and Vel_Error is not Zero then Change_Velocity is High"
  "if Position_Error is not Zero and Vel_Error is High then Change_Velocity is High"
];

% ruleset_accel = [
%   "if Position_Error is Zero and Vel_Error is Zero then Accel_Out is Zero"
%   "if Position_Error is Moderate or Vel_Error is Moderate then Accel_Out is Moderate"
%   "if Position_Error is Zero and Vel_Error is High then Accel_Out is Moderate"
%   "if Position_Error is High and Vel_Error is Zero then Accel_Out is Moderate"
%   "if Position_Error is High and Vel_Error is High then Accel_Out is High"
% ];

% ruleset_velocity = [
%   "if Position_Error is Zero and Vel_Error is not High then Change_Velocity is Zero"   
%   "if Position_Error is Zero and Vel_Error is High then Change_Velocity is Moderate"    
%   "if Position_Error is not Zero and Vel_Error is Zero then Change_Velocity is Moderate"   
%   "if Position_Error is not Zero and Vel_Error is Moderate then Change_Velocity is Moderate"    
%   "if Position_Error is not Zero and Vel_Error is High then Change_Velocity is High"
% ];

ruleset_accel = [
  "if Position_Error is not High and Vel_Error is Zero then Accel_Out is Zero"   
  "if Vel_Error is Moderate then Accel_Out is Moderate"   
  "if Position_Error is not High and Vel_Error is High then Accel_Out is Moderate"   
  "if Position_Error is High and Vel_Error is not High then Accel_Out is Moderate"   
  "if Position_Error is High and Vel_Error is High then Accel_Out is High"
];

acc_fis = addRule(acc_fis, ruleset_velocity);
acc_fis = addRule(acc_fis, ruleset_accel);

clear ruleset_velocity ruleset_accel


% Custom gensurf where we only look at the normalized inputs of the error
%   terms.

num_x = 30; num_y = 30; num_squares = num_x * num_y;

% Create a custom grid with specified ranges
[input1, input2] = meshgrid(        ...
      linspace(0, PE_width, num_x), ...
      linspace(0, VE_width, num_y)  ...
);

in_1_flat  = reshape(input1, num_squares, 1);
in_2_flat  = reshape(input2, num_squares, 1);

% graphing 2d surfaces ----------------------------
Z       = evalfis( acc_fis, [in_1_flat, in_2_flat] );
Z1_flat = reshape( Z(:,1), num_y, num_x );
Z2_flat = reshape( Z(:,2), num_y, num_x );

f1 = figure; f1.Position(1:2) = sp(3,1); hold on
  surf(input1, input2, Z1_flat);
  xlabel 'Position Error [m]';
  ylabel 'Velocity Error [m/s]';
  zlabel 'Change of Velocity [m/s]';
  title  'Change of Velocity Surface Graph'

f1 = figure; f1.Position(1:2) = sp(3,2); hold on
  surf(input1, input2, Z2_flat);
  xlabel 'Position Error [m]';
  ylabel 'Velocity Error [m/s]';
  zlabel 'Accel_Out [m/s^2]';
  title  'Acceleration-Output Surface Graph'


clear Z Z1_flat Z2_flat in_1_flat in_2_flat input1 input2


pause_tour(tour_each_subsys, "fuzzy speed controller is done.")





% Model Preparations -----------------------------------------------------
window = @(t, a,b)   (t >= a) - (t >= b);  % window function


% time properties ---------------------------------
dt        =  1/5;   % seconds
max_time  =  30;    % seconds

time = 0:dt:max_time; % seconds

g = 9.81; % m/s^2


% other cars --------------------------------------
speed_2       = 30;   % m/s
position_2    = 100;  % m


% uncomment to set the constant-speed case
% x_car_2 =  (position_2 + time.*speed_2);


% details for setting up the decreasing-speed case
w   = -10;  % controls acceleration
t2  =  10;  % seconds; after this time, speed decreases to zero


% uncomment to set the decreasing-speed-till-stop case
x_car_2 =  (position_2 + time.*speed_2) .* window(time,  0, t2)     ...
         + parabola(                                                ...
             time, w,                                               ...
             t2,   speed_2, position_2) .* window(time, t2, t2 - w) ...
         + 550.*window(time, t2-w, max_time);  % m

x_car_2(end) = 550;     % needs to be adjusted when 'w' and 't2' change

x_car_2(1)  = position_2;


% This initial distance has to be set correctly or else the controller will
%   go wild.
distance_from_front = rain_intensity * evalfis(bd_fis, speed_2); % m


% This will be filled in as we go through the control loop
set_point     = zeros(1, length(time));
vel_set_pt    = zeros(1, length(time));
accel_set_pt  = zeros(1, length(time));


% state space model (SSM) ---------------------------------------------
% the inputs of the SSM include:
% U_n = [ a_out;
%         v_out;]
%   a_out and v_out are also the outputs of the fuzzy speed controller


% state vector  -----------------------
% previous X_n      |   next X_n
% X_n   = [ a_n;    |   X_np1 = [ a_np1;
%           v_n;    |             v_np1;
%           x_n ];  |             x_np1 ];

% initial state
X_0 = [ 0;
        speed_2;
        position_2 - distance_from_front  ];

X_n = X_0;

A = [ 0         0   0;
      0         1   0;
      0.5*dt^2  dt  1  ];

B = [ 1         0;
      dt        1;
      0.5*dt^2  dt];

% output is 3. this means it is controllable to any state X_n that we want.
rank_svm = rank(ctrb(A,B));


% generalized state variable model format:
% X_np1 = A*X_n + B*U_n


graph_array_1 = zeros(length(time), 4);
graph_array_2 = zeros(length(time), 3);


% control/animation loop -------------------------------------------------


disable_animation = true;

% animation 
if disable_animation == false

  fig       = figure; hold on;
  myWriter  = VideoWriter('Car_no_rain - 10 mps case.mp4', 'MPEG-4');
  myWriter.FrameRate = 10;
  open(myWriter);

  % car parameters
  carWidth = 20; carHeight = 6;
  
  wheel_radius  =  2;
  xoffset       = [3, 16];
  yoffset       =  3;
  
  roof = carWidth./3;
  roofColor = [0.2, 0.2, 0.2]; % RGB for darker gray
  
  % define multiple cars
  numCars           =  2;
  car_positions     = [0, 40];
  colors            = ['b', 'r'];
  
  % Creating car
  cars    = cell(1, numCars);
  wheels  = cell(1, numCars, 2);
  roofs   = cell(1, numCars);
  
  theta   = linspace(0, 2*pi, 50); % Angle values for circles
  wheelX  = wheel_radius * cos(theta);
  wheelY  = wheel_radius * sin(theta);
  
  semicircle  = linspace(0, pi, 50);
  roofx       = roof * cos(semicircle);
  roofy       = roof * sin(semicircle);
  
  % rainy parameters
  rain              = .6;           % 60% of normal performance
  % max_acceleration  = 15 * rain;    % reducing the max acceleration by 60%
  % max_deceleration  = 15.5 * rain;  % reducing the max deceleration by 60%
  % safe_distance     = 10/rain;      % increasing safe distance in meters
  % set_speed         = 15 * rain;    % reducing speed in m/s 
  
  % Rain
  rain_density  = 400; % # of rain drops
  rain_drops    = gobjects(1, rain_density);
  
  % axis limit
  x_range = [0 1500];
  y_range = [0 50];
  
  for i = 1:rain_density
    x = rand() * diff(x_range) + x_range(1);
    y = rand() * diff(y_range) + y_range(1);
    rain_drops(i) = plot(x,y, '.', 'Color', [0.5, 0.5, 1], 'MarkerSize', 10);
  end
  
  
  for i = 1:numCars
    % Create wheels and roofs
    cars{i}       = rectangle(                                            ...
            'Position', [car_positions(i), yoffset, carWidth, carHeight], ...
            'FaceColor', colors(i)                                        ...
    );
    wheels{i, 1}  = fill(car_positions(i) + xoffset(1) + wheelX, yoffset + wheelY, 'black');
    wheels{i, 2}  = fill(car_positions(i) + xoffset(2) + wheelX, yoffset + wheelY, 'black');
  
    roofs{i}      = fill(car_positions(i) + carWidth/2 + roofx, yoffset + carHeight + roofy, roofColor);
  end

end


% Control loop --------------------------------------------------------
for n = 1:length(time)
  n;

  % This is where rain_intensity is changed at a certain sample number
  if n >= 20 && is_it_raining
    rain_intensity = rain_intensity + 0.05;
    rain_intensity = limiter(rain_intensity, 1,2);
  end
  

  % Step 1: given a previous velocity input, the distance_from_front is
  %     evaluated.
  if mod(n, bd_update_rate) == 0
    vel_in = X_n(2);
    distance_from_front = rain_intensity * evalfis(bd_fis, vel_in);
  end


  % Step 2: we update our set point
  set_point(n+1)      = x_car_2(n)        - distance_from_front;
  vel_set_pt(n+1)     = (set_point(n+1)   - set_point(n))/dt;
  accel_set_pt(n+1)   = (vel_set_pt(n+1)  - vel_set_pt(n))/dt;
  
  set_pt_n = [ accel_set_pt(n+1);
                 vel_set_pt(n+1);
                  set_point(n+1)  ];
  

  % Step 3: Get the summing block output, which is the error term.
  error_n = set_pt_n - X_n;


  % Step 4: Given the two error terms, we evaluate the improvement terms.
  %     Here, we only consider the magnitude of each of the error terms.
  %                               pos_error   vel_error
  output = evalfis( acc_fis, abs([error_n(3), error_n(2)]) );
  
  
  % If x_error is (+) the vel_output should be (-), while accel_output ...
  %               (-)                          (+)
  %   ... should point opposite of vel_output.
  out_accel   =  sign(error_n(3)) * output(1);
  out_vel     = -sign(error_n(3)) * output(2);
  

  % Step 5: initialize the input_n vector in the SSM format
  input_n = [ out_accel;
              out_vel    ];


  % Step 6: then, we evaluate the state vector X_n+1
  X_np1   = A*X_n + B*input_n;
  

  % place the vals onto the handy graph_array
  graph_array_1(n, 5) = X_np1(1);     % acceleration
  graph_array_1(n, 4) = X_np1(2);     % velocity
  graph_array_1(n, 2) = X_np1(3);     % position
  graph_array_1(n, 3) = set_pt_n(3);  
  graph_array_1(n, 1) = x_car_2(n);

  graph_array_2(n, 1) = rain_intensity;
  graph_array_2(n, 2) = distance_from_front;
  graph_array_2(n, 3) = x_car_2(n) - X_np1(3);
  

  % optional: print tenth result ------------------
  if n == 10 && false
    graph_array_1
  end
  

  % Step 7: update animation ---------------------------------------------
  if disable_animation == false

    car_positions(1) = X_np1(3) - carWidth;
    car_positions(2) = x_car_2(n);


    for i = 1:numCars
      set(  cars{i},                                                      ...
            'Position', [car_positions(i), yoffset, carWidth, carHeight]  ...
      );
  
      % Update wheels position for car i
      set(  wheels{i, 1},                                     ...
            'XData', car_positions(i) + xoffset(1) + wheelX,  ...
            'YData', yoffset + wheelY                         ...
      );
      set(                                                ...
        wheels{i, 2},                                     ...
        'XData', car_positions(i) + xoffset(2) + wheelX,  ...
        'YData', yoffset + wheelY                         ...
      );
  
      % Update roof position for car i
      set(  roofs{i}, ...
            'XData', car_positions(i) + carWidth/2 + roofx, ...
            'YData', yoffset + carHeight + roofy ...
      );
    end

    if rain_intensity > 1
      for i = 1:rain_density
        x = get(rain_drops(i), 'XData');
        y = get(rain_drops(i), 'YData') - 1; % creates rain falling

        if y < y_range(1) % reset rain
          y = y_range(2);
          x = rand() * diff(x_range) + x_range(1);
        end

        set(rain_drops(i), 'XData', x, 'YData', y);
      end
    end
    
    % adjusting axis to follow the two cars
    xlim([ min(car_positions) - 30,  max(car_positions) + 30 ]);
    ylim(y_range)

    % capture frames
    frame = getframe(gcf);
    writeVideo(myWriter, frame);
  
    pause(.05); % Pause for a short duration to create animation effect

  end


  % Step 7: update the previous state variable X_n with the new one for
  %   the time when sample time 'dt' passes. This is the 1/z operator!
  X_n = X_np1;

  % When the file plays, the last sample throws a "fuzzy.internal ...
  %   ... .utility.throwWarning" warning. This is because the window
  %   function isn't extending all the way to the very last sample. But it
  %   can effectively be ignored.
end

if disable_animation == false
  close(myWriter);
end


pause_tour(tour_each_subsys, "simulation is done.")



% post-animation calculations --------------------------------------------

x_error = graph_array_1(:, 2) - graph_array_1(:, 3);

% sum of the square of error
MSE_total = sum(x_error.^2)

MSE_t = zeros(length(time), 1);
for k = 1:(length(time) - 1)
  MSE_t(k+1) = MSE_t(k) + x_error(k)^2;
end



% Graphs Galore ----------------------------------------------------------


f1 = figure; f1.Position(1:2) = sp(2,1); hold on
  yyaxis left, stairs(time, graph_array_2(:,1), LineWidth=2), ylim padded
    title 'Objective Parameters', xlabel 'Time [seconds]'
    ylabel('Rain Intensity')
  yyaxis right, stairs(time, graph_array_2(:,2:3), LineWidth=2), ylim padded
    ylabel('Distance [m]'), grid on
    legend({'rain intensity',                         ...
            'braking distance set pt.',                           ...
            'actual'}, Location="northwest")


% graphing positions of all objects and position set point across time
graphing_array = graph_array_1(:, 1:3);

f1 = figure; f1.Position(1:2) = sp(4,1); hold on
subplot(7,1, [1 5]), stairs(time, graphing_array, LineWidth=2)
  title 'Positions [m] of each Object over Time'
  legend({'front car Position',                         ...
          'ACC Car Position',                           ...
          'positional set point'}, Location="northwest")
  ylabel('Position [m]'), grid on


% graphing velocity and acceleration across time
accel = graph_array_1(:, 5);
vel   = graph_array_1(:, 4);

subplot(7,1, [6 7]), hold on, grid on, xlabel('Time [seconds]')
  yyaxis right, stairs(time, accel, LineWidth=2), ylim padded
    ylabel 'Acceleration [m/s^2]'
  yyaxis left, stairs(time, vel, LineWidth=2), ylim padded
    ylabel 'Velocity [m/s]'


% graphing error and MSE across time
f1 = figure; f1.Position(1:2) = sp(4,2); hold on, grid on
  title 'Error across Time'
  yyaxis left, stairs(time, x_error, LineWidth=2)
    ylabel 'Error [m]'
  yyaxis right, stairs(time, MSE_t, LineWidth=2)
    ylabel 'MSE'

  legend({'error[n]',                                     ...
          ['MSE[n]; total error = ' num2str(MSE_total)]}, ...
          Location="west")
  
  xlabel 'Time [seconds]'




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User Defined Functions (UDFs) ------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pause_tour(bool, text)
  if bool
    fprintf("\n")
    fprintf(text)
    fprintf("\npress any key in the command window to continue --------------------------")
    pause
  end
  fprintf("\n\n\n")
  close all
end


function out = parabola(time, w, t2, speed_2, position_2)
  % details for setting up the increasing- or decreasing-speed case
  
  % parabola function: x = a(t-h)^2 + k
  a   = speed_2/(2*w);
  h   = t2 - w;
  k   = speed_2*(t2 - w/2) + position_2;
  
  out = a*(time - h).^2 + k;  % m
end