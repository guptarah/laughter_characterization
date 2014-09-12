function NonHomPP = poisson_analysis()

% create poisson vector
win_size = 10;
files01 = dir('../../data/to_use_set_01/*.turn'); 

all_data = [];
for file = files01'
	cur_file=strcat('../../data/to_use_set_01/',file.name);
	cur_data = load(cur_file);
	max_index = floor(size(cur_data,1)/win_size)*win_size;

	all_data = [all_data; cur_data(1:max_index,:)];
end

counselor_laughs = (all_data(:,1) == 0).*all_data(:,3);
client_laughs = (all_data(:,1) == 1).*all_data(:,3);

Config.PhysioFs=1;
Config.CellWidth=30; 
Config.Alpha=0.05;
Config.NumOfIterations=100; 
Config.RateFunctionName='rectwinsum';
Config.TimeUnit=10;

addpath('../PoissonProcess/');

% first doing for client laughs
[HomPP] = EstimateAndGenerateHomogeneousPoissonProcess(client_laughs',Config); % homogeneous estimation
ChangePoints = find(counselor_laughs);
Config.ChangePoints=ChangePoints;
[Config.InitialParameters] = InitializePoissonRateFunctionParameters(Config.RateFunctionName,length(Config.ChangePoints),HomPP.EstimatedParameters);
[NonHomPP] = EstimateAndGeneratePoissonProcess2(client_laughs',Config);

disp('for client laughs');
HomPP
NonHomPP

% doing it for counselor laughs
disp('for counselor laughs');
[HomPP] = EstimateAndGenerateHomogeneousPoissonProcess(counselor_laughs',Config); % homogeneous estimation
HomPP
ChangePoints = find(client_laughs);
Config.ChangePoints=ChangePoints; 
[Config.InitialParameters] = InitializePoissonRateFunctionParameters(Config.RateFunctionName,length(Config.ChangePoints),HomPP.EstimatedParameters);
[NonHomPP] = EstimateAndGeneratePoissonProcess2(counselor_laughs',Config);
NonHomPP 
