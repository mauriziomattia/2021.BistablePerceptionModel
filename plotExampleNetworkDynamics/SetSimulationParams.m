function Net = SetSimulationParams(Net,Params)

X1 = 1; %Population indexing
X2 = 2;
Y1 = 3;
Y2 = 4;

TauIntegration= Params(1); %Time-constant, integration layer
TauSwitch = Params(2); %Time-constant, decision layer
WE_self = Params(3);
WI_top = Params(4);
WE_ff = Params(5);
WI_fb = Params(6);
WI_AS = Params(7);
OffsetSwitch = Params(8);

Net.CParam.W(Y2,X1) = -WI_AS; % Effective connectivity matrix
Net.CParam.W(X2,Y1) = -WI_AS;
Net.CParam.W(X2,X1) = WE_ff-WI_AS;
Net.CParam.W(Y2,Y1) = WE_ff-WI_AS;
Net.CParam.W(X2,X2) = WE_self;
Net.CParam.W(Y2,Y2) = WE_self;
Net.CParam.W(X1,X2) = -WI_fb;
Net.CParam.W(Y1,Y2) = -WI_fb;
Net.CParam.W(X2,Y2) = -WI_top;
Net.CParam.W(Y2,X2) = -WI_top;


% Offsets and characteristic time constants
Net.SPParam.Theta(X1) = 0; % The contrast...
Net.SPParam.Theta(X2) = OffsetSwitch;
Net.SPParam.Tau0(X1) = TauIntegration;
Net.SPParam.Tau0(X2) = TauSwitch;
Net.SPParam.Theta(Y1) = 0; % The contrast...
Net.SPParam.Theta(Y2) = OffsetSwitch;
Net.SPParam.Tau0(Y1) = TauIntegration;
Net.SPParam.Tau0(Y2) = TauSwitch;

Net.n(X1) = 0;
Net.n(X2) = 0;
Net.n(Y1) = 0;
Net.n(Y2) = 0;

Net.Adapt.r0 = 0; % No need for this (adaptation on top inhibition)
Net.Adapt.TauSTD = 0.5;

return