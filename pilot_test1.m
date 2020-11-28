clear VCSW BLSW BOSWT TUSW1 TUSW2 TUSW3 TUSW4 SW1 SW2 SW3 SW4

VCSW = VertCompletionObj;
VCSW.Pws = 2e4;
VCSW.Twf = 45;
VCSW.WellPIJ = 20/100;

VCSW1 = VertCompletionObj;
VCSW1.Pws = 2e4;
VCSW1.Twf = 45;
VCSW1.WellPIJ = 25/100;

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

BOSWT = BOObj;
BOSWT.Tbosc = 851 - 460;
BOSWT.Tbgsc = 196.47 - 460;
BOSWT.gSG = 0.64;
BOSWT.wSG = 1.00;
BOSWT.oSG = GE_oil(2);
BOSWT.GOR = 130;
BOSWT.WC = 0;
BOSWT.TuneViscosity(BLSW.bmrmu, BLSW.bmrT*1.8 + 32)

TUSW1 = TubingObj;
TUSW1.TVD = 100;
TUSW1.tTgprof = - 0.025;
TUSW1.tln = 100;
TUSW1.tdi = 0.1053;
TUSW1.trg = 1.524e-5;
TUSW1.tuc = 11.349/1e3;
TUSW1.tin = 0;

TUSW2 = TubingObj;
TUSW2.TVD = 500;
TUSW2.tTgprof = - 0.025;
TUSW2.tln = 500;
TUSW2.tdi = 0.1053;
TUSW2.trg = 1.524e-5;
TUSW2.tuc = 11.349/1e3;
TUSW2.tin = 0;

TUSW3 = TubingObj;
TUSW3.TVD = 300;
TUSW3.tTgprof = - 0.025;
TUSW3.tln = 500;
TUSW3.tdi = 0.1053;
TUSW3.trg = 1.524e-5;
TUSW3.tuc = 11.349/1e3;
TUSW3.tin = 0;

TUSW4 = TubingObj;
TUSW4.TVD = 200;
TUSW4.tTgprof = - 0.025;
TUSW4.tln = 500;
TUSW4.tdi = 0.1053;
TUSW4.trg = 1.524e-5;
TUSW4.tuc = 11.349/1e3;
TUSW4.tin = 0;

%PPSW2 = ESPObj;
%PPSW2.CurveHQ = [H2 Q2];
%PPSW2.BEP = [60*21.1 2*151.7];
%PPSW2.N =  3600;
%PPSW2.Nc = 3600;

SW1 = SingleWellObj;
SW1.mrT = BLSW.bmrT;
SW1.mrmu = BLSW.bmrmu;
SW1.BOmodel = BOSWT;
SW1.items = {VCSW TUSW1};
SW1.WC = 0;

SW2 = SingleWellObj;
SW2.mrT = BLSW.bmrT;
SW2.mrmu = BLSW.bmrmu;
SW2.BOmodel = BOSWT;
SW2.items = {VCSW TUSW2};
SW2.WC =0;

SW3 = SingleWellObj;
SW3.mrT = BLSW.bmrT;
SW3.mrmu = BLSW.bmrmu;
SW3.BOmodel = BOSWT;
SW3.items = {VCSW TUSW3};
SW3.WC = 0;

SW4 = SingleWellObj;
SW4.mrT = BLSW.bmrT;
SW4.mrmu = BLSW.bmrmu;
SW4.BOmodel = BOSWT;
SW4.items = {VCSW TUSW4};
SW4.WC = 0;

FL1 = FlowlineObj;
FL1.tTinf = 15;
FL1.tln = 100;
FL1.tdi = .0254*10;
FL1.trg = 1.524e-5;
FL1.tuc = 11.349/1e3;
FL1.CalculationType = 'forward';

FL2 = FlowlineObj;
FL2.tTinf = 15;
FL2.tln = 200;
FL2.tdi = .0254*10;
FL2.trg = 1.524e-5;
FL2.tuc = 11.349/1e3;
FL2.CalculationType = 'forward';

FL3 = FlowlineObj;
FL3.tTinf = 15;
FL3.tln = 200;
FL3.tdi = .0254*10;
FL3.trg = 1.524e-5;
FL3.tuc = 11.349/1e3;
FL3.CalculationType = 'forward';

FL4 = FlowlineObj;
FL4.tTinf = 15;
FL4.tln = 100;
FL4.tdi = .0254*10;
FL4.trg = 1.524e-5;
FL4.tuc = 11.349/1e3;
FL4.CalculationType = 'forward';

FL5 = FlowlineObj;
FL5.tTinf = 15;
FL5.tln = 100;
FL5.tdi = .0254*10;
FL5.trg = 1.524e-5;
FL5.tuc = 11.349/1e3;
FL5.CalculationType = 'forward';

FL6 = FlowlineObj;
FL6.tTinf = 15;
FL6.tln = 300;
FL6.tdi = .0254*10;
FL6.trg = 1.524e-5;
FL6.tuc = 11.349/1e3;
FL6.CalculationType = 'forward';

B1 = {SW1 FL1};
B2 = {SW2 FL2};
B3 = {SW3 FL3};
B4 = {SW4 FL4};

N1 = NodeObj;
N1.BOmodel = BOSWT;
N1.Branches = {B1 B2};

bconnectivity = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
nconnectivity = [0 1 0; 0 0 0; 0 1 0];

%SW.items = {VCSW TUSW1};

