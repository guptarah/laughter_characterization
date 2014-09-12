rng('shuffle','multFibonacci');

% configuration parameters
Config.PhysioFs=32;
Config.CellWidth=30;
Config.Alpha=0.05;
Config.NumOfIterations=100;
Config.RateFunctionName='expsum';
Config.TimeUnit=5;

% load signal that we want to model with a PP
load('SCROccurencesExample.mat','SCROccurences');

% add random change points in rate function
load('ChangePoints','ChangePoints');

% Homogeneous PP
[HomPP] = EstimateAndGenerateHomogeneousPoissonProcess(SCROccurences,Config);

% Non-homogeneous PP - Parametric rate function
Config.ChangePoints=ChangePoints;
[Config.InitialParameters] = InitializePoissonRateFunctionParameters(Config.RateFunctionName,length(Config.ChangePoints),HomPP.EstimatedParameters);
[NonHomPP] = EstimateAndGeneratePoissonProcess2(SCROccurences,Config);

% Non-homogeneous PP - Piecewise linear non-parametric rate function
[PieceLinPP] = EstimateAndGeneratePoissonProcessPiecewiseLinear(SCROccurences,Config);