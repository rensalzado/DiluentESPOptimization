clear VCSW1 VCSW2 VCSW3 BLSW1 BLSW2 BLSW3 BOSW1 BOSW2 BOSW3...
      TUSW1 TUSW2 TUSW3 TUSW4 TUSW5 TUSW6 ...
      PPSW ...
      SW1 SW2 SW3 N1 N2 N3
%% Vertical completion
VCSW1 = VertCompletionObj;
VCSW1.Pws = 2e4;
VCSW1.Twf = 45;
VCSW1.WellPIJ = 30/100;

VCSW2 = VertCompletionObj;
VCSW2.Pws = 2e4;
VCSW2.Twf = 55;
VCSW2.WellPIJ = 28/100;

VCSW3 = VertCompletionObj;
VCSW3.Pws = 2e4;
VCSW3.Twf = 50;
VCSW3.WellPIJ = 20/100;

%% Injection point
BLSW1 = InjectionObj;
BLSW1.Qo = 0;
BLSW1.Qd = 500;
BLSW1.Qw = 0;
BLSW1.Tb = 30;
BLSW1.mrmu = [vmu_oil(2, :); vmu_dil(3, :)]; 
BLSW1.mrT = [Tsm_oil(2, :); Tsm_dil(3, :)];
BLSW1.drho = GE_dil(3)*1e3;
BLSW1.orho = GE_oil(2)*1e3;
BLSW1.BlendingReference

BLSW2 = InjectionObj;
BLSW2.Qo = 0;
BLSW2.Qd = 500;
BLSW2.Qw = 0;
BLSW2.Tb = 30;
BLSW2.mrmu = [vmu_oil(3, :); vmu_dil(3, :)]; 
BLSW2.mrT = [Tsm_oil(3, :); Tsm_dil(3, :)];
BLSW2.drho = GE_dil(3)*1e3;
BLSW2.orho = GE_oil(3)*1e3;
BLSW2.BlendingReference

BLSW3 = InjectionObj;
BLSW3.Qo = 0;
BLSW3.Qd = 500;
BLSW3.Qw = 0;
BLSW3.Tb = 30;
BLSW3.mrmu = [vmu_oil(3, :); vmu_dil(3, :)]; 
BLSW3.mrT = [Tsm_oil(3, :); Tsm_dil(3, :)];
BLSW3.drho = GE_dil(3)*1e3;
BLSW3.orho = GE_oil(3)*1e3;
BLSW3.BlendingReference

BLSW01 = InjectionObj;
BLSW01.Qo = 1;
BLSW01.Qd = 1;
BLSW01.Qw = 0;
BLSW01.Tb = 40;
BLSW01.mrmu = [vmu_oil(2, :); vmu_oil(2, :)]; 
BLSW01.mrT = [Tsm_oil(2, :); Tsm_oil(2, :)];
BLSW01.drho = GE_oil(2)*1e3;
BLSW01.orho = GE_oil(2)*1e3;
BLSW01.BlendingReference

BLSW02 = InjectionObj;
BLSW02.Qo = 1;
BLSW02.Qd = 1;
BLSW02.Qw = 0;
BLSW02.Tb = 40;
BLSW02.mrmu = [vmu_oil(3, :); vmu_oil(3, :)]; 
BLSW02.mrT = [Tsm_oil(3, :); Tsm_oil(3, :)];
BLSW02.drho = GE_oil(3)*1e3;
BLSW02.orho = GE_oil(3)*1e3;
BLSW02.BlendingReference

BLSW03 = InjectionObj;
BLSW03.Qo = 1;
BLSW03.Qd = 1;
BLSW03.Qw = 0;
BLSW03.Tb = 40;
BLSW03.mrmu = [vmu_oil(2, :); vmu_oil(2, :)]; 
BLSW03.mrT = [Tsm_oil(2, :); Tsm_oil(2, :)];
BLSW03.drho = GE_oil(2)*1e3;
BLSW03.orho = GE_oil(2)*1e3;
BLSW03.BlendingReference

%% Black oil model
BOSW1 = BOObj;
BOSW1.Tbosc = 851 - 460;
BOSW1.Tbgsc = 196.47 - 460;
BOSW1.gSG = 0.64;
BOSW1.wSG = 1.00;
BOSW1.oSG = GE_oil(2);
BOSW1.GOR = 130;
BOSW1.WC = 0;
BOSW1.TuneViscosity(BLSW01.bmrmu, BLSW01.bmrT*1.8 + 32)

BOSW2 = BOObj;
BOSW2.Tbosc = 851 - 460;
BOSW2.Tbgsc = 196.47 - 460;
BOSW2.gSG = 0.64;
BOSW2.wSG = 1.00;
BOSW2.oSG = GE_oil(3);
BOSW2.GOR = 130;
BOSW2.WC = 0;
BOSW2.TuneViscosity(BLSW02.bmrmu, BLSW02.bmrT*1.8 + 32)

BOSW3 = BOObj;
BOSW3.Tbosc = 851 - 460;
BOSW3.Tbgsc = 196.47 - 460;
BOSW3.gSG = 0.64;
BOSW3.wSG = 1.00;
BOSW3.oSG = GE_oil(3);
BOSW3.GOR = 130;
BOSW3.WC = 0;
BOSW3.TuneViscosity(BLSW03.bmrmu, BLSW03.bmrT*1.8 + 32)

%% Tubing 
% SW1
% section 1
TUSW1 = TubingObj;
TUSW1.TVD = 100;
TUSW1.tTgprof = -0.025;
TUSW1.tln = 100;
TUSW1.tsl = 5;
TUSW1.tdi = 0.1053;
TUSW1.trg = 1.524e-5;
TUSW1.tuc = 11.349/1e3;
TUSW1.tin = 0;
% section 2
TUSW2 = TubingObj;
TUSW2.TVD = 300;
TUSW2.tTgprof = -0.025;
TUSW2.tln = 300;
TUSW2.tsl = 15;
TUSW2.tdi = 0.1053;
TUSW2.trg = 1.524e-5;
TUSW2.tuc = 11.349/1e3;
TUSW2.tin = 0;

% SW2
% section 1
TUSW3 = TubingObj;
TUSW3.TVD = 100;
TUSW3.tTgprof = -0.025;
TUSW3.tln = 100;
TUSW3.tsl = 5;
TUSW3.tdi = 0.1053;
TUSW3.trg = 1.524e-5;
TUSW3.tuc = 11.349/1e3;
TUSW3.tin = 0;

% section 2
TUSW4 = TubingObj;
TUSW4.TVD = 300;
TUSW4.tTgprof = -0.025;
TUSW4.tln = 300;
TUSW4.tsl = 15;
TUSW4.tdi = 0.1053;
TUSW4.trg = 1.524e-5;
TUSW4.tuc = 11.349/1e3;
TUSW4.tin = 0;

% SW3
% section 1
TUSW5 = TubingObj;
TUSW5.TVD = 100;
TUSW5.tTgprof = -0.025;
TUSW5.tln = 100;
TUSW5.tsl = 5;
TUSW5.tdi = 0.1053;
TUSW5.trg = 1.524e-5;
TUSW5.tuc = 11.349/1e3;
TUSW5.tin = 0;

% section 2
TUSW6 = TubingObj;
TUSW6.TVD = 300;
TUSW6.tTgprof = -0.025;
TUSW6.tln = 300;
TUSW6.tsl = 15;
TUSW6.tdi = 0.1053;
TUSW6.trg = 1.524e-5;
TUSW6.tuc = 11.349/1e3;
TUSW6.tin = 0;

%% ESP 2
PPSW = ESPObj;
PPSW.CurveHQ = [H2 Q2];
PPSW.BEP = [30*21.1 3*151.7];
PPSW.N =  3600;
PPSW.Nc = 3600;

%% SW configuration
% SW1
SW1 = SingleWellObj;
SW1.mrT = BLSW01.bmrT;
SW1.mrmu = BLSW01.bmrmu;

SW1.BOmodel = BOSW1;
SW1.items = {VCSW1 TUSW1 BLSW1 PPSW TUSW2};

% SW2
SW2 = SingleWellObj;
SW2.mrT = BLSW02.bmrT;
SW2.mrmu = BLSW02.bmrmu;

SW2.BOmodel = BOSW2;
SW2.items = {VCSW2 TUSW3 BLSW2 PPSW TUSW4};

% SW3
SW3 = SingleWellObj;
SW3.mrT = BLSW03.bmrT;
SW3.mrmu = BLSW03.bmrmu;

SW3.BOmodel = BOSW3;
SW3.items = {VCSW3 TUSW5 BLSW3 PPSW TUSW6};

%% FL configuration

FL1 = FlowlineObj;
FL1.tTinf = 15;
FL1.tln = 300;
FL1.tsl = 15;
FL1.tdi = .1053;
FL1.trg = 1.524e-5;
FL1.tuc = 11.349/1e3;
FL1.CalculationType = 'forward';

FL2 = FlowlineObj;
FL2.tTinf = 15;
FL2.tln = 200;
FL2.tsl = 10;
FL2.tdi = .1053;
FL2.trg = 1.524e-5;
FL2.tuc = 11.349/1e3;
FL2.CalculationType = 'forward';

FL3 = FlowlineObj;
FL3.tTinf = 15;
FL3.tln = 50;
FL3.tsl = 5;
FL3.tdi = .254;
FL3.trg = 1.524e-5;
FL3.tuc = 11.349/1e3;
FL3.CalculationType = 'forward';

FL4 = FlowlineObj;
FL4.tTinf = 15;
FL4.tln = 300;
FL4.tsl = 20;
FL4.tdi = .254;
FL4.trg = 1.524e-5;
FL4.tuc = 11.349/1e3;
FL4.CalculationType = 'forward';

FL5 = FlowlineObj;
FL5.tTinf = 15;
FL5.tln = 300;
FL5.tsl = 20;
FL5.tdi = .254;
FL5.trg = 1.524e-5;
FL5.tuc = 11.349/1e3;
FL5.CalculationType = 'forward';

%% Node configuration

B1 = {SW1 FL1};
B2 = {SW2 FL2};
B3 = {SW3 FL3};

N1 = NodeObj;
N1.P = 7000;
N1.Branches = {B1 B2};

N2 = NodeObj;
N2.P = 7000;
N2.Branches = {B3};

NB1 = {N1 FL4};
NB2 = {N2 FL5};

N3 = NodeObj;
N3.P = 6500;
N3.Nodes = {NB1 NB2};