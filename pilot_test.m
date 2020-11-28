%% Vertical completion
VCSW = VertCompletionObj;
VCSW.Pws = 2e4;
VCSW.Twf = 45;
VCSW.WellPIJ = 30/100;

%% Injection point
BLSW = InjectionObj;
BLSW.Qo = 500;
BLSW.Qd = 0;
BLSW.Qw = 0;
BLSW.Tb = 30;
BLSW.mrmu = [vmu_oil(2, :); vmu_dil(3, :)]; 
BLSW.mrT = [Tsm_oil(2, :); Tsm_dil(3, :)];
BLSW.drho = GE_oil(2)*1e3;
BLSW.orho = GE_dil(3)*1e3;
BLSW.BlendingReference

%% Black oil model
BOSW = BOObj;
BOSW.Tbosc = 851 - 460;
BOSW.Tbgsc = 196.47 - 460;

BOSW.gSG = 0.64;
BOSW.wSG = 1.00;
BOSW.oSG = GE_oil(2);
BOSW.GOR = 130;
BOSW.WC = 0;
BOSW.TuneViscosity(BLSW.bmrmu, BLSW.bmrT*1.8 + 32)

%% Tubing standard
TUSW = TubingObj;
TUSW.TVD = 500;
TUSW.tTgprof = 0.1;
TUSW.tln = 500;
TUSW.tdi = 0.1053;
TUSW.trg = 1.524e-5;
TUSW.tuc = 11.349/1e3;
TUSW.tin = 0;

%% Tubing first section
TUSW1 = TubingObj;
TUSW1.TVD = 100;
TUSW1.tTgprof = 0.2;
TUSW1.tln = 100;
TUSW1.tdi = 0.1053;
TUSW1.trg = 1.524e-5;
TUSW1.tuc = 11.349/1e3;
TUSW1.tin = 0;

%% Tubing second section
TUSW2 = TubingObj;
TUSW2.TVD = 500;
TUSW2.tTgprof = 0.2;
TUSW2.tln = 500;
TUSW2.tdi = 0.1053;
TUSW2.trg = 1.524e-5;
TUSW2.tuc = 11.349/1e3;
TUSW2.tin = 0;

%% ESP 1
PPSW = ESPObj;
PPSW.CurveHQ = [H1 Q1];
PPSW.BEP = [1228.2 454.3];
PPSW.N =  3600;
PPSW.Nc = 2900;

%% ESP 2
PPSW2 = ESPObj;
PPSW2.CurveHQ = [H2 Q2];
PPSW2.BEP = [60*21.1 2*151.7];
PPSW2.N =  3600;
PPSW2.Nc = 3600;


%% SW configuration 1
SW = SingleWellObj;
SW.mrT = BLSW.bmrT;
SW.mrmu = BLSW.bmrmu;

SW.BOmodel = BOSW;
SW.items = {VCSW TUSW};

%% SW configuration 2
SW = SingleWellObj;
SW.mrT = BLSW.bmrT;
SW.mrmu = BLSW.bmrmu;

SW.BOmodel = BOSW;
SW.items = {VCSW BLSW TUSW};

%% SW configuration 3
SW = SingleWellObj;
SW.mrT = BLSW.bmrT;
SW.mrmu = BLSW.bmrmu;

SW.BOmodel = BOSW;
SW.items = {VCSW TUSW1 BLSW PPSW2 TUSW2};


