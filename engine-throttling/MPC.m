nlobj = nlmpc(1,1,1);

nlobj.Ts = 0.01;
nlobj.PredictionHorizon = 80;
nlobj.ControlHorizon = 10;

nlobj.Model.StateFcn = "engineStateFcn";
nlobj.Model.OutputFcn = "engineOutputFcn";
nlobj.Model.IsContinuousTime = true;

mu = 0;       
sigma = 0.025;    
n = 1000;     
noise = mu + sigma * randn(n, 1);


nlobj.MV.Min = 0;
nlobj.MV.Max = 0.55 + noise;
Max = nlobj.MV.Max;

nlobj.MV.RateMin = -0.4;
nlobj.MV.RateMax = 0.4;

nlobj.Weights.OutputVariables = 20;
nlobj.Weights.ManipulatedVariables = 0;
nlobj.Weights.ManipulatedVariablesRate = 0.001;