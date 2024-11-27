%% NOS CHILL MAIN
% File to perform the main iterative calculations for the refrigeration
% optimization effort - minimize tchill with constraints on power input
% this is a tool to design the heat transfer and thermodynamic behaviour of
% a NOS chilling system. Due consideration must be given to general system
% design, and these inputs can be iterated within this program to find the
% optimal system

% refrigeration cycles occur at constant pressures 
% 1-2: isentropic compression from Pevap to Pcond
% 2-3: isobaric heat rejection in condenser to saturated liquid at Pcond
% 3-4: isenthalpic throttle from Pcond to Pevap
% 4-1: isobaric heat absorption in evaporator to saturated vapour at Pevap

%% TO ASK
% how to deal with flow of two phase refrigerant
    % particularly pressure drop, calculating heat transfer coefficiencts
    % how to deal with this across the expansion valve
% how to deal with convection from two phase substance

% he says there are empirical correlations, or EES might have a thing
% split it up into convection from the liquid and convection from vapour


%% Inputs, Ambient and Target Conditions
% Tamb = ambient temperature, starting NOS temperature
% Ttarget = desired NOS temperature, saturation conditions
% Vnos = volume of the tank, assume filled with NOS
% mnos = mass of nitrous in the tank

% Modules
    % CoolProp OR RefProp (NIST based, purhcase required)

%% Optimization Goals and Ouputs
% t = time of refrigeration (minimize)
% W = power draw (minimize)

% Variables
    % mdotRef = refrigerant mass flow rate
    % Dpipe = pipe diameter for refrigerant
    % Tevap1 = first evaporator temperature for stage 1
    % Tevap2 is defined by being below the final desired NOS temperature
    % Tcond1 = first condenser temperature for stage 1
    % Tcond2 is defined by being above the ambient temperature

% Plots


%% Design Requirements
% Heat transfer mechanism
    % produce the most effective heat transfer mechanism (high rate)
    % look at cooling channel design (areas, flow rates, etc)
    % look at insulation design in the outer shell
    % look at transient heat transfer, flowing fluids, state change, etc.
% Thermodynamic efficiency
    % design the most efficienct thermodyanmic cycle for the working fluid
    % dual stage cycle may be needed for low temperature
    % select fluids with appropriate bp temperature for operating range
% Material selection
    % select appropriate tube materials, sizes, etc
% NOS thermodynamic behaviour
    % transient behaviour of chilling nitrous oxide
% Power system bahviour, compatibility, etc.
    % power requirements, power supply design, etc


%% Heat Transfer Analysis
% have a working fluid refrigerant, which is an design variable
% heat absorption occurs at the evaporator temperature
    % here, phase change of the refrigerant from liquid to vapour occurs
% heat rejection occurs at the condenser temperature
    % here, phase change of the refrigerant from vapour to liquid occurs
% mass flow rate should match the cycle with the rate of heat transfer
% in this case, working temperatures are optimization variables
% mass flow rate is a result of the heat transfer calculation

% EVAPORATION - time stepped
% Start condition: refrigerant evaporator pressure
    % refrigerant at saturated liquid condition, know T12
% Process: latent heat absorption and evaporation at constant pressure
    % adjust for temperature by calling properties in time stepping
% End condition: refrigerant evaporator pressure
    % refrigerant at superheated vapour condition, T changed
% Calculated: total heat absorption by refrigerant
% Back calculated: rate of heat transfer through all systems
    % have a total heat transferred in one pass
    % need to calculate change in NOS temperature with transient heat flux
    % NOS final temperature is iteration RESULT (iterative solution)
    % assume a final temperature of NOS in one pass
    % calculate heat transfer to working fluid in one pass by iterating
    % several timesteps backward from this temperature
    % remember the temperature adjustment of teh refrigerant with each time
    % step, accounted for with the mass flow rate
    % if it does not match the energy for the evaporation process, repeat
% Result: final temperature of NOS after one pass of refrigerant


% COMPRESSOR - single power calculation
% optimization variables Pevap and Pcond
% can directly calculate work by the compressor with refrigerant states
% account for losses here


% CONDENSATION - time stepped
% Start condition: refrigerant condenser temperature
    % refrigerant at liquid-vap or vapour condition depending on cycle
% Process: latent heat rejection and condensation at constant temperature
% End condition: refrigerant condenser temperature
    % refrigerant at saturated? state
% Calculated: total heat rejection by refrigerant
% Back calculated: system requirements (length, flow rate, etc. to achieve)
    % have total heat trasnferred in one pass
    % calculate transient heat flux to environment
    % calculate required length, etc for this
    % no converging iteration required, only time stepping
% Result: system design requirements after one pass of refrigerant


% EXPANSION VALVE - single calculation for Cv
% optimization variables Pevap and Pcond
% mass flow rate known
% calculate required valve Cv to achieve this pressure drop
% assume losses here deviating from (will affect the 4-1 curve slope)
% make a constant loss assumption so states are all known, simplify it


% OVERALL - loop stepped
% based on optimization variables Pevap, Pcond, process cycle is known
% W for compressor is directly calculated
% repeat loops of evaporation and condensation until temp is reached
% try to minimize time to cool
% take some overall optimization function that combines W and t in a
% meaningful way to minimize and give the best design


% DEVIATONS
% flow losses create a pressure drop in the system
% creates deviations in heat transfer properties of the working fluid
% needs to be augmented with the time stepping
% flow losses will be created in the valve and compressor
% this will further mess everything up


%% Project Outcomes
% Optimization model for working pressures
% Optimization tests (manually) for heat exchanger sizing, design
% Optimization tests (manually) for working fluid, stages
% Material selection for heat exchanger, seals, etc.
% Pressure sims for tanks, piping, etc.
% P&ID of heat exchanger system
% Cost analysis of system
% Power analysis of system
% CAD of full design including tanks and plumbing
% Full heat transfer analysis of final system and insulation


%% FUNCTIONS
%% HEAT EXCHANGER FLOW
% Flow modelling for the design of a in-tank heat exchanger
% Goal - maximize surface area
% Concept - bulkhead mounted and inserted into nitrous tank before fill
%% Refrigerant Convection Integrated - UPDATE
% calculates total resistance of the refrigerant as a lengthwise integral
% convergent solution required similar to other convection cases
function [Rref,dPf,ho] = RefrigerantConvection(Ti,si,Tpipei,Lpipe,dt,dl,Dpipei,mdot,Pevap,tol,D,q)
    % inlet refrigerant states
    T = Ti;
    s = si;
    Pi = py.CoolProp.CoolProp.PropsSI('P', 'T', T, 'S', s, 'R134a');

    % initialize guess value and differental length
    Ts = Tpipei;
    dh = q/mdot;
    hi = py.CoolProp.CoolProp.PropsSI('H', 'T', Ti, 'S', si, 'R134a');
    ho = hi + dh;
    havg = (hi+ho)/2;
    

    x = py.CoolProp.CoolProp.PropsSI('Q', 'P', Pi, 'HMASS', havg, 'R134a');


    % compute superficial mass flux G
    G = mdot/(pi/4*(Dpipei^2));

    % compute Reynolds, Relo and Revo
    muf = py.CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 0, 'R134a');
    mug = py.CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 1, 'R134a');
    kf = py.CoolProp.CoolProp.PropsSI('L', 'T', T, 'Q', 0, 'R134a');           % Thermal conductivity [W/m.K]
    cpf = py.CoolProp.CoolProp.PropsSI('C', 'T', T, 'Q', 0, 'R134a');          % Specific heat capacity [J/kg.K]
    rhog = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 1, 'R134a');
    rhof = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 0, 'R134a');
    k = py.CoolProp.CoolProp.PropsSI('L', 'T', T, 'Q', 0, 'R134a');
    sgma = py.CoolProp.CoolProp.PropsSI('I', 'T', T, 'Q', 0, 'R134a');

    % Compute dimensionless numbers
    Relo = G*Dpipei*(1-x)/muf;
    Revo = G*Dpipei*x/mug;
    Prf = muf*cpf/kf;
    De = (Relo+Revo*(mug/muf)*((rhof/rhog)^0.5))*((Dpipei/D)^0.5);

    % compute Lockhart-Martinelli parameter
    LMparam = (((1-x)/x)^0.9)*((rhog/rhof)^0.5)*((muf/mug)^0.1);

    % compute enthalpy
    hfg = py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 1, 'R134a') - py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 0, 'R134a');

    % compute relevant heat transfer areas for convection to pipe
    Aref = Lpipe*pi*Dpipei;

    % compute boiling number Boi
    Boi = q/(Aref*G*hfg);

    % Compute Nusselt number
    Nu = 7850*(De^0.43)*(Prf^(-5.055))*((Boi*(10^4))^0.125)*(LMparam^(-0.036));

    % compute local h
    h = k*Nu/Dpipei;

    % compute q correction
    Rref = 1/(h*Aref);

    % pressure drop
    % Friedel
    if Relo < 2000
        flo = 16/Relo;
    else
        flo = 0.079*(Relo^(-0.25));
    end
    if Revo < 2000
        fvo = 16/Revo;
    else
        fvo = 0.079*(Revo^(-0.25));
    end
    g = 9.81;
    E = ((1-x)^2)+(x^2)*(rhof/rhog)*(fvo/flo);
    F = (x^0.76)*((1-x)^0.24);
    H = ((rhof/rhog)^0.91)*((mug/muf)^0.19)*((1-(mug/muf))^0.7);
    rhohom = (x/rhog + (1-x)/rhof)^(-1);
    Fr = (G^2)/(g*Dpipei*(rhohom^2));
    We = ((G^2)*Dpipei)/(sgma*rhohom);

    multiplier2 = E + 3.23*(F*H/((Fr^0.043)*(We^0.035)));
    dPf = flo*(2*dl/Dpipei)*((G^2)/(g*rhof))*multiplier2;

    P = Pi-dPf;

    % update state
    To = py.CoolProp.CoolProp.PropsSI('T', 'H', ho, 'P', P, 'R134a');
    T = (Ti+To)/2;
    P = (Pi+P)/2;

    % recompute
    x = py.CoolProp.CoolProp.PropsSI('Q', 'P', P, 'H', havg, 'R134a');

    % compute superficial mass flux G
    G = mdot/(pi/4*(Dpipei^2));

    % compute Reynolds, Relo and Revo
    muf = py.CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 0, 'R134a');
    mug = py.CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', 1, 'R134a');
    kf = py.CoolProp.CoolProp.PropsSI('L', 'T', T, 'Q', 0, 'R134a');           % Thermal conductivity [W/m.K]
    cpf = py.CoolProp.CoolProp.PropsSI('C', 'T', T, 'Q', 0, 'R134a');          % Specific heat capacity [J/kg.K]
    rhog = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 1, 'R134a');
    rhof = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 0, 'R134a');
    k = py.CoolProp.CoolProp.PropsSI('L', 'T', T, 'Q', 0, 'R134a');
    sgma = py.CoolProp.CoolProp.PropsSI('I', 'T', T, 'Q', 0, 'R134a');

    % Compute dimensionless numbers
    Relo = G*Dpipei*(1-x)/muf;
    Revo = G*Dpipei*x/mug;
    Prf = muf*cpf/kf;
    De = (Relo+Revo*(mug/muf)*((rhof/rhog)^0.5))*((Dpipei/D)^0.5);

    % compute Lockhart-Martinelli parameter
    LMparam = (((1-x)/x)^0.9)*((rhog/rhof)^0.5)*((muf/mug)^0.1);

    % compute enthalpy
    hfg = py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 1, 'R134a') - py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 0, 'R134a');

    % compute relevant heat transfer areas for convection to pipe
    Aref = Lpipe*pi*Dpipei;

    % compute boiling number Boi
    Boi = q/(Aref*G*hfg);

    % Compute Nusselt number
    Nu = 7850*(De^0.43)*(Prf^(-5.055))*((Boi*(10^4))^0.125)*(LMparam^(-0.036));

    % compute local h
    h = k*Nu/Dpipei;

    % compute q correction
    Rref = 1/(h*Aref);

    % pressure drop
    % Friedel
    if Relo < 2000
        flo = 16/Relo;
    else
        flo = 0.079*(Relo^(-0.25));
    end
    if Revo < 2000
        fvo = 16/Revo;
    else
        fvo = 0.079*(Revo^(-0.25));
    end
    g = 9.81;
    E = ((1-x)^2)+(x^2)*(rhof/rhog)*(fvo/flo);
    F = (x^0.76)*((1-x)^0.24);
    H = ((rhof/rhog)^0.91)*((mug/muf)^0.19)*((1-(mug/muf))^0.7);
    rhohom = (x/rhog + (1-x)/rhof)^(-1);
    Fr = (G^2)/(g*Dpipei*(rhohom^2));
    We = ((G^2)*Dpipei)/(sgma*rhohom);

    multiplier2 = E + 3.23*(F*H/((Fr^0.043)*(We^0.035)));
    dPf = flo*(2*dl/Dpipei)*((G^2)/(g*rhof))*multiplier2;

    P = Pi-dPf;
end

%% NOS CHILL THERMODYNAMICS
% File including functions related to instantaneous thermodynamic relations
% Tnos = NOS temperature
% Vnos = NOS tank volume
% Lnos = NOS tank length
% mnos = NOS mass

%% NOS Level
% liquid level of nitrous at given state
function [Lliq,Lvap] = NOSLevel(Tnos,Vnos,mnos,Lnos,Dtanki)
    rhovap = py.CoolProp.CoolProp.PropsSI('D','T',Tnos,'Q',1,'CO2');
    x = py.CoolProp.CoolProp.PropsSI('Q','T',Tnos,'D',mnos/Vnos,'CO2');
    mvap = x*mnos;
    Vvap = mvap/rhovap;
    Lvap = Vvap/(pi/4*(Dtanki^2));
    Lliq = Lnos -Lvap;
end

%% Exchanger Level
% length of heat exchanger immersed in liquid and vapour
function [Lpliq,Lpvap] = PipeLevel(Lpipe,Hpipe,Lliq)
    % compute the pipe length to height ratio
    LtoH = Lpipe/Hpipe;

    % compute the length immersed in liquid and vapour
    Lpliq = Lliq*LtoH;
    Lpvap = Lpipe - Lpliq;
end

%% Refrigeration Cycle States - UPDATE WITH NON-IDEALITIES
% computes states of the refrigeration cycle given input pressures
function [s,h,T] = RefrigerationCycle(Pevap,Pcond)
    % state 2 - pre-throttle
    T2 = py.CoolProp.CoolProp.PropsSI('T','P',Pcond,'Q',0,'R134a');
    h2 = py.CoolProp.CoolProp.PropsSI('H','P',Pcond,'Q',0,'R134a');
    s2 = py.CoolProp.CoolProp.PropsSI('S','P',Pcond,'Q',0,'R134a');

    % state 3 - pre-evaporator
    h3 = h2;
    T3 = py.CoolProp.CoolProp.PropsSI('T','P',Pevap,'H',h3,'R134a');
    s3 = py.CoolProp.CoolProp.PropsSI('S','P',Pevap,'H',h3,'R134a');  

    % state 4 - post-evaporator
    T4 = py.CoolProp.CoolProp.PropsSI('T','P',Pevap,'Q',1,'R134a');
    h4 = py.CoolProp.CoolProp.PropsSI('H','P',Pevap,'Q',1,'R134a');
    s4 = py.CoolProp.CoolProp.PropsSI('S','P',Pevap,'Q',1,'R134a');    

    % state 1 - post-compressor
    s1 = s4;
    T1 = py.CoolProp.CoolProp.PropsSI('T','P',Pcond,'S',s1,'R134a');
    h1 = py.CoolProp.CoolProp.PropsSI('H','P',Pcond,'S',s1,'R134a');

    % states
    s = [s1,s2,s3,s4];
    h = [h1,h2,h3,h4];
    T = [T1,T2,T3,T4];
end

%% NOS CHILL HEAT TRANSFER
% File including all functions related to instantaneous heat transfer rates
% assume perfectly insulated top and bottom of tank

%% NOS Convection Pipe - CHECK
% convection with two phase NOS at instantaneous temperature, proportion
% assume saturation conditions maintained, P and T change simultaneously
function [Rnosliq,Rnosvap] = NOSConvectionPipe(Tnos,TpipeO,Lpliq,Lpvap,DpipeO,Hpipe,n,p,D,inc,Lpipe)
    % compute film temperature guess
    Tf = (Tnos + TpipeO)/2;
    HtoL = Hpipe/Lpipe;
    % compute Ra and Pr for free convection, use D
    % liquid properties
    rhof = py.CoolProp.CoolProp.PropsSI('D', 'T', Tf, 'Q', 0, 'CO2');         % Density [kg/m^3]
    mu = py.CoolProp.CoolProp.PropsSI('V', 'T', Tf, 'Q', 0, 'CO2');          % Dynamic viscosity [Pa.s]
    kliq = py.CoolProp.CoolProp.PropsSI('L', 'T', Tf, 'Q', 0, 'CO2');           % Thermal conductivity [W/m.K]
    cp = py.CoolProp.CoolProp.PropsSI('C', 'T', Tf, 'Q', 0, 'CO2');          % Specific heat capacity [J/kg.K]
    
    % Calculating derived properties
    nuf = mu / rhof;                                        % Kinematic viscosity [m^2/s]
    beta = py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', 'T', Tf, 'Q', 0, 'CO2');  % Volumetric thermal expansion coefficient [1/K]
    alpha = kliq / (rhof * cp);                               % Thermal diffusivity [m^2/s]
    g = 9.81;

    % Rayleight and Prandtl numbers, and H/D parameter
    RaliqH = g*beta*(Tnos-TpipeO)*((Lpliq*HtoL)^3)/(nuf*alpha);
    RaliqL = g*beta*(Tnos-TpipeO)*((Lpliq)^3)/(nuf*alpha);
    RaliqD = g*beta*(Tnos-TpipeO)*(D^3)/(nuf*alpha);
    Prliq = mu*cp/kliq;
    HD = (n+1)*DpipeO/D + n*p/D;

    % vapour properties
    rhog = py.CoolProp.CoolProp.PropsSI('D', 'T', Tf, 'Q', 1, 'CO2');         % Density [kg/m^3]
    mu = py.CoolProp.CoolProp.PropsSI('V', 'T', Tf, 'Q', 1, 'CO2');          % Dynamic viscosity [Pa.s]
    kvap = py.CoolProp.CoolProp.PropsSI('L', 'T', Tf, 'Q', 1, 'CO2');           % Thermal conductivity [W/m.K]
    cp = py.CoolProp.CoolProp.PropsSI('C', 'T', Tf, 'Q', 1, 'CO2');          % Specific heat capacity [J/kg.K]
    
    % Calculating derived properties
    nug = mu / rhog;                                        % Kinematic viscosity [m^2/s]
    beta = py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', 'T', Tf, 'Q', 1, 'CO2');  % Volumetric thermal expansion coefficient [1/K]
    alpha = kvap / (rhog * cp);                               % Thermal diffusivity [m^2/s]
    g = 9.81;

    % Prandtl number, Ja, B, and hfg for vapour
    Prvap = mu*cp/kvap;
    hfg = py.CoolProp.CoolProp.PropsSI('H', 'T', Tnos, 'Q', 1, 'N2O') - py.CoolProp.CoolProp.PropsSI('H', 'T', Tnos, 'Q', 0, 'N2O');
    Ja = cp*(Tnos-TpipeO)/hfg;
    hfgprime = hfg*(1+(0.683-0.228/Prvap)*Ja);
    B = (rhof-rhog)/rhof*(cp*(Tnos-TpipeO)/hfgprime)*((tan(inc))^2)/Prvap;
    
    % define Z table
    Blist = [0,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
    dDlist = [0.1,0.2,0.3,0.4,0.5];
    Zlist = [1,1,1,1,1;
            1,1.001,1.001,1.003,1.004;
            1,1.001,1.003,1.005,1.008;
            1.001,1.002,1.004,1.008,1.012;
            1.001,1.003,1.006,1.010,1.017;
            1.001,1.003,1.007,1.013,1.021;
            1.002,1.007,1.015,1.028,1.047;
            1.002,1.010,1.023,1.045,1.078;
            1.003,1.013,1.032,1.063,1.117;
            1.004,1.017,1.042,1.085,1.167;
            1.005,1.021,1.052,1.109,1.235;
            1.006,1.025,1.062,1.136,1.329;
            1.007,1.029,1.073,1.167,1.451;
            1.008,1.033,1.085,1.201,1.557;
            1.008,1.037,1.098,1.238,1.601];
    Z = interp2(dDlist,Blist,Zlist,DpipeO/D,B);

    % compute Nu using empirical correlation
    if (RaliqH > 1.49E5 && RaliqH < 1.43E6)
        Nuliq = ((RaliqH^0.2925)*(HD^(-0.1515)));
    elseif (RaliqL > 5.01E10 && RaliqL < 7E11)
        Nuliq = ((RaliqL^0.3071)*(HD^(-0.1097)));
    %elseif (RaliqD > 1.29E6 && RaliqD < 3.04E8)
    else
        Nuliq = ((RaliqD^0.2039)*(HD^(-0.1798)));
    end

    % compute h from Nu for liq and vap
    hliq = Nuliq*kliq/DpipeO;
    hvap = (g*(((rhof*(rhof-rhog))*hfgprime*(kliq^3)*cos(inc))/(3*DpipeO*nug*(Tnos-TpipeO)*Z)))^0.25;

    % compute relevant heat transfer areas for convection to pipe
    Apliq = Lpliq*pi*DpipeO;
    Apvap = Lpvap*pi*DpipeO;

    % compute liquid and vapour resistance coefficients of convection
    Rnosliq = 1/(hliq*Apliq);
    Rnosvap = 1/(hvap*Apvap);
end

%% Piping Conduction - GOOD
% conduction with instantaneous piping resistance
function Rpipe = PipingConduction(Lpipe,DpipeO,Dpipei,kpipe)
    % compute required pipe geometries
    Rpipe2 = DpipeO/2;
    Rpipe1 = Dpipei/2;

    % compute 1D pipe resistance
    Rpipe = log(Rpipe2/Rpipe1)/(2*pi*Lpipe*kpipe);
end

%% NOS Convection Tank - UPDATE (IF IN FREE CONVECTION REGIME)
% convection with two phase NOS at instantaneous temperature, proportion
% assume saturation conditions maintained, P and T change simultaneously
function [Rnosliq,Rnosvap] = NOSConvectionTank(Tnos,Ttanki,Lpliq,Lpvap,Dtanki)
    % compute film temperature guess
    Tf = (Tnos + Ttanki)/2;

    % compute Ra and Pr for free convection, use D
    % liquid properties
    rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Tf, 'Q', 0, 'CO2');         % Density [kg/m^3]
    mu = py.CoolProp.CoolProp.PropsSI('V', 'T', Tf, 'Q', 0, 'CO2');          % Dynamic viscosity [Pa.s]
    kliq = py.CoolProp.CoolProp.PropsSI('L', 'T', Tf, 'Q', 0, 'CO2');           % Thermal conductivity [W/m.K]
    cp = py.CoolProp.CoolProp.PropsSI('C', 'T', Tf, 'Q', 0, 'CO2');          % Specific heat capacity [J/kg.K]
    
    % Calculating derived properties
    nu = mu / rho;                                        % Kinematic viscosity [m^2/s]
    beta = py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', 'T', Tf, 'Q', 0, 'CO2');  % Volumetric thermal expansion coefficient [1/K]
    alpha = kliq / (rho * cp);                               % Thermal diffusivity [m^2/s]
    g = 9.81;

    % Rayleight and Prandtl numbers
    Raliq = g*beta*(Ttanki-Tnos)*(Dtanki^3)/(nu*alpha);
    Prliq = mu*cp/kliq;

    % vapour properties
    rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Tf, 'Q', 1, 'CO2');         % Density [kg/m^3]
    mu = py.CoolProp.CoolProp.PropsSI('V', 'T', Tf, 'Q', 1, 'CO2');          % Dynamic viscosity [Pa.s]
    kvap = py.CoolProp.CoolProp.PropsSI('L', 'T', Tf, 'Q', 1, 'CO2');           % Thermal conductivity [W/m.K]
    cp = py.CoolProp.CoolProp.PropsSI('C', 'T', Tf, 'Q', 1, 'CO2');          % Specific heat capacity [J/kg.K]
    
    % Calculating derived properties
    nu = mu / rho;                                        % Kinematic viscosity [m^2/s]
    beta = py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', 'T', Tf, 'Q', 1, 'CO2');  % Volumetric thermal expansion coefficient [1/K]
    alpha = kvap / (rho * cp);                               % Thermal diffusivity [m^2/s]
    g = 9.81;

    % Rayleight and Prandtl numbers
    Ravap = g*beta*(Ttanki-Tnos)*(Dtanki^3)/(nu*alpha);
    Prvap = mu*cp/kvap;

    % compute Nu using empirical correlation (vertical pipe free) UPDATE
    Nuliq = (0.825 + (0.387*(Raliq^(1/6)))/(1+(0.492/Prliq)^(9/16))^(8/27))^2;
    Nuvap = (0.825 + (0.387*(Ravap^(1/6)))/(1+(0.492/Prvap)^(9/16))^(8/27))^2;

    % compute h from Nu for liq and vap
    hliq = Nuliq*kliq/Dtanki;
    hvap = Nuvap*kvap/Dtanki;

    % compute relevant heat transfer areas for convection to pipe
    Apliq = Lpliq*pi*Dtanki;
    Apvap = Lpvap*pi*Dtanki;

    % compute liquid and vapour resistance coefficients of convection
    Rnosliq = 1/(hliq*Apliq);
    Rnosvap = 1/(hvap*Apvap);
end

%% Tank Conduction - GOOD
% conduction with instantaneous tank temperature gradient
function Rtank = TankConduction(Lnos,DtankO,Dtanki,ktank)
    % compute required geometries
    Rtank2 = DtankO/2;
    Rtank1 = Dtanki/2;

    % compute 1D cylinder resistance
    Rtank = log(Rtank2/Rtank1)/(2*pi*Lnos*ktank);
end

%% Contact Resistance
% given R''
function Rcont = ContactResistance(Lnos,DtankO,Rpp)
    % compute required geometry
    Acont = pi*DtankO*Lnos;

    % compute Rcont
    Rcont = Rpp/Acont;
end

%% Insulation Conduction - GOOD
% 1D cylinder condution
function Rins = InsulationConduction(Lins,Dins,DtankO,kins)
    % compute required geometries
    Rins2 = Dins/2;
    Rins1 = DtankO/2;

    % compute 1D cylinder resistance
    Rins = log(Rins2/Rins1)/(2*pi*Lins*kins);
end

%% Ambient Convection - GOOD
% assume stationary environment fluid
function Ramb = AmbientConvection(Tamb,TinsO,Lins,Dins)
    % compute film temperature guess
    Tf = (Tamb + TinsO)/2;
    P = 101325;

    % compute Ra and Pr for free convection, use D
    % liquid properties
    rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Tf, 'P', P, 'Air');         % Density [kg/m^3]
    mu = py.CoolProp.CoolProp.PropsSI('V', 'T', Tf, 'P', P, 'Air');          % Dynamic viscosity [Pa.s]
    k = py.CoolProp.CoolProp.PropsSI('L', 'T', Tf, 'P', P, 'Air');           % Thermal conductivity [W/m.K]
    cp = py.CoolProp.CoolProp.PropsSI('C', 'T', Tf, 'P', P, 'Air');          % Specific heat capacity [J/kg.K]
    
    % Calculating derived properties
    nu = mu / rho;                                        % Kinematic viscosity [m^2/s]
    beta = py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', 'T', Tf, 'P', P, 'Air');  % Volumetric thermal expansion coefficient [1/K]
    alpha = k / (rho * cp);                               % Thermal diffusivity [m^2/s]
    g = 9.81;

    % Rayleight and Prandtl numbers
    Ra = g*beta*(Tamb-TinsO)*(Lins^3)/(nu*alpha);
    Pr = mu*cp/k;

    % compute Nu using empirical correlation (vertical wall free)
    Nu = (0.825 + (0.387*(Ra^(1/6)))/(1+(0.492/Pr)^(9/16))^(8/27))^2;

    % compute h from Nu
    hair = Nu*k/Lins;

    % compute resistance coefficient of convection
    Ramb = 1/(hair*pi*Dins*Lins);
end

%% Heat From NOS - UPDATE
% solve for Ttest = Tnos to get q out of NOS
function [Ttest,dPfc,hoc] = HeatOut(q,Ti,To,si,so,Tnos,Vnos,mnos,Lnos,Dtanki,Lpipe,Hpipe,Dpipei,DpipeO,kpipe,mdot,tol,dt,dl,Pevap,D,n,p,inc)
    % compute levels
    [Lliq,Lvap] = NOSLevel(Tnos,Vnos,mnos,Lnos,Dtanki);
    [Lpliq,Lpvap] = PipeLevel(Lpipe,Hpipe,Lliq);
    Tavg = (Ti+To)/2;

    % individual temperatures with qguess, for Ttest (NOS)
    % refrigerant to inner wall
    %Rref = RefrigerantConvection(Ti,To,si,so,Lpipe,Dpipei,mdot,tol);
    toln = tol + 1;
    a = Ti;
    b = Tnos;
    iterations = 0;
    while abs(toln) > tol
        c = (a + b)/2;
        Ts = c;
        if iterations > 20
            tol = tol+0.1;
        end
        [Rrefc,dPfc,hoc] = RefrigerantConvection(Ti,si,c,Lpipe,dt,dl,Dpipei,mdot,Pevap,tol,D,q)
        [Rrefa,dPfa,hoa] = RefrigerantConvection(Ti,si,a,Lpipe,dt,dl,Dpipei,mdot,Pevap,tol,D,q)
        if ((Ti+q*Rrefc)-c)*((Ti+q*Rrefa)-a) < 0
            b = c
        else
            a = c
        end
        toln = (Ti+q*Rrefc)-c
        tol
        iterations = iterations + 1;
    end
    Tpipei = Ts;

    % inner wall to outer wall
    Rpipe = PipingConduction(Lpipe,DpipeO,Dpipei,kpipe);
    TpipeO = q*Rpipe + Tpipei;

    if TpipeO > Tnos
        Ttest = Tnos+100;
    else
        % outer wall to NOS
        [Rliqpipe,Rvappipe] = NOSConvectionPipe(Tnos,TpipeO,Lpliq,Lpvap,DpipeO,Hpipe,n,p,D,inc,Lpipe);
        Rnos = 1/(1/Rliqpipe+1/Rvappipe);
        Ttest = q*Rnos + TpipeO;
    end
end

%% Heat From Ambient - UPDATE
% neglect radiation
function Ttest = HeatIn(q,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)
    % compute levels
    [Lliq,Lvap] = NOSLevel(Tnos,Vnos,mnos,Lnos,Dtanki);

    % individual temperatures with qguess, for Ttest (NOS)
    % ambient to outer insulation
    toln = tol + 1;
    a = Tnos;
    b = Tamb;
    while abs(toln) > tol
        c = (a + b)/2;
        Rambc = AmbientConvection(Tamb,c,Lins,Dins);
        Ramba = AmbientConvection(Tamb,a,Lins,Dins);
        if ((Tamb-q*Rambc)-c)*((Tamb-q*Ramba)-a) < 0
            b = c;
        else
            a = c;
        end
        toln = (Tamb-q*Rambc)-c;
    end
    TinsO = c;

    % outer insulation to inner insulation
    Rins = InsulationConduction(Lins,Dins,DtankO,kins);
    Tinsi = TinsO - q*Rins;

    % contact resistance
    Rcont = ContactResistance(Lnos,DtankO,Rpp);
    TtankO = Tinsi - q*Rcont;

    % outer tank to inner tank
    Rtank = TankConduction(Lnos,DtankO,Dtanki,ktank);
    Ttanki = TtankO - q*Rtank;

    if Ttanki < Tnos
        Ttest = Tnos-100
    else
        % inner tank to NOS
        [Rliqtank,Rvaptank] = NOSConvectionTank(Tnos,Ttanki,Lliq,Lvap,Dtanki);
        Rnos = 1/(1/Rliqtank+1/Rvaptank);
    
        Ttest = Ttanki - q*Rnos;
    end
end

%% HEAT TRANSFER - UPDATE
% transient heat transfer analysis for NOS cooling
function [Tnos,x,qin,qout,Pnos,dPfn,hon] = HeatTransfer(dt,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,DpipeO,Ti,To,si,so,Lpipe,Hpipe,Dpipei,kpipe,mdot,qnsin,qnsout,tol,tolq,dl,Pevap,D,n,p,inc)
    % guess q, calculate all temperatures for in and out
    % compute heat addition rate from ambient to NOS
    toln = tolq + 1;
    qns0 = qnsin;
    qns1 = qns0-0.001;
    iterations = 0;
    while abs(toln) > tolq
        fns0 = HeatIn(qns0,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        fns1 = HeatIn(qns1,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        qns = qns1 - (fns1*(qns1-qns0)/(fns1-fns0))
        % if toln == -Inf
        %     qns0 = 0.01
        %     qns1 = 0.02;
        %     disp('reset')
        % else
        if iterations > 20 || (fns1 == -100 && fns0 == -100)
            qns0 = 0.001
            qns1 = 0.002
            tol = tol+1
            iterations = 0;
        %     qns = qnssaved
        %     tol = HeatIn(qns,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos + 0.1
        %     toln = tol
        else
            toln = HeatIn(qns,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
            qns0 = qns1
            qns1 = qns
        %     if fns0 ~= -100 && toln ~= -100
        %         qnssaved = qns0
        %         tolsaved = HeatIn(qnssaved,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        %     elseif fns1 ~= -100 && toln ~= -100
        %         qnssaved = qns1
        %         tolsaved = HeatIn(qnssaved,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        %     end
            iterations = iterations + 1;
        end
    end
    qin = qns;

    % compute heat extraction rate from liquid and vapour to refrigerant
    toln = tolq + 1;
    qns0 = qnsout;
    qns1 = qns0+0.1;

    dh0 = qns0/mdot
    while dh0 > 38280
        qns0 = qns0/2;
    end
    dh1 = qns1/mdot
    while dh1 > 38280
        qns1 = qns1/2;
    end

    iterations = 0;
    while abs(toln) > tolq
        [Tns0,dPf0,ho0] = HeatOut(qns0,Ti,To,si,so,Tnos,Vnos,mnos,Lnos,Dtanki,Lpipe,Hpipe,Dpipei,DpipeO,kpipe,mdot,tol,dt,dl,Pevap,D,n,p,inc);
        fns0 = Tns0 - Tnos;
        [Tns1,dPf1,ho1] = HeatOut(qns1,Ti,To,si,so,Tnos,Vnos,mnos,Lnos,Dtanki,Lpipe,Hpipe,Dpipei,DpipeO,kpipe,mdot,tol,dt,dl,Pevap,D,n,p,inc);
        fns1 = Tns1 - Tnos;
        qns = qns1 - (fns1*(qns1-qns0)/(fns1-fns0))

        dh = qns/mdot;
        if dh > 38280
            qns0 = qns0/2;
            qns1 = qns1/2;
            iterations = 0;
        else
            if iterations > 20
                tol = tol+1;
            end
            dh
            qns
            [Tns,dPfn,hon] = HeatOut(qns,Ti,To,si,so,Tnos,Vnos,mnos,Lnos,Dtanki,Lpipe,Hpipe,Dpipei,DpipeO,kpipe,mdot,tol,dt,dl,Pevap,D,n,p,inc);
            toln = Tns - Tnos
            qns0 = qns1
            qns1 = qns
            iterations = iterations + 1;
        end
    end
    qout = qns
    toln

    % compute heat extracted in dt = change in internal energy
    dQ = qout*dt - qin*dt;

    % compute new internal energy of NOS
    u0 = py.CoolProp.CoolProp.PropsSI('U','T',Tnos,'D',mnos/Vnos,'CO2');
    u1 = (u0*mnos - dQ)/mnos;

    % update NOS state at fixed density, new internal energy
    % assume equilibrium is re-established quickly
    Tnos = py.CoolProp.CoolProp.PropsSI('T','U',u1,'D',mnos/Vnos,'CO2');
    Pnos = py.CoolProp.CoolProp.PropsSI('P','U',u1,'D',mnos/Vnos,'CO2');
    x = py.CoolProp.CoolProp.PropsSI('Q','U',u1,'D',mnos/Vnos,'CO2');
end

function [Tnos,x,qin,Pnos] = IdleHeating(dt,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,DpipeO,Ti,To,si,so,Lpipe,Hpipe,Dpipei,kpipe,mdot,qnsin,qnsout,tol,tolq,dl,Pevap,D,n,p,inc)
    % guess q, calculate all temperatures for in and out
    % compute heat addition rate from ambient to NOS
    toln = tolq + 1;
    qns0 = qnsin;
    qns1 = qns0-0.001;
    iterations = 0;
    while abs(toln) > tolq
        fns0 = HeatIn(qns0,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        fns1 = HeatIn(qns1,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        qns = qns1 - (fns1*(qns1-qns0)/(fns1-fns0))
        % if toln == -Inf
        %     qns0 = 0.01
        %     qns1 = 0.02;
        %     disp('reset')
        % else
        if iterations > 20 || (fns1 == -100 && fns0 == -100)
            qns0 = 0.001
            qns1 = 0.002
            tol = tol+1
            iterations = 0;
        %     qns = qnssaved
        %     tol = HeatIn(qns,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos + 0.1
        %     toln = tol
        else
            toln = HeatIn(qns,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
            qns0 = qns1
            qns1 = qns
        %     if fns0 ~= -100 && toln ~= -100
        %         qnssaved = qns0
        %         tolsaved = HeatIn(qnssaved,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        %     elseif fns1 ~= -100 && toln ~= -100
        %         qnssaved = qns1
        %         tolsaved = HeatIn(qnssaved,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,tol)-Tnos
        %     end
            iterations = iterations + 1;
        end
    end
    qin = qns;

    % compute heat extracted in dt = change in internal energy
    dQ = -qin*dt;

    % compute new internal energy of NOS
    u0 = py.CoolProp.CoolProp.PropsSI('U','T',Tnos,'D',mnos/Vnos,'CO2');
    u1 = (u0*mnos - dQ)/mnos;

    % update NOS state at fixed density, new internal energy
    % assume equilibrium is re-established quickly
    Tnos = py.CoolProp.CoolProp.PropsSI('T','U',u1,'D',mnos/Vnos,'CO2');
    Pnos = py.CoolProp.CoolProp.PropsSI('P','U',u1,'D',mnos/Vnos,'CO2');
    x = py.CoolProp.CoolProp.PropsSI('Q','U',u1,'D',mnos/Vnos,'CO2');
end

%% POWER
% function to compute power and total energy consumed with pressure losses
function [W,WP,dE,dEP] = Power(Pcond,Pi,dPf,mdot,ho,eta,dt)
    % compute conditions at evaporator outlet
    Po = Pi - dPf;
    so = py.CoolProp.CoolProp.PropsSI('S','P',Po,'H',ho,'R134a');
    ho

    % isentropic efficiency on compressor - assume very low (L+V)
    sconds = so;
    hconds = py.CoolProp.CoolProp.PropsSI('H','P',Pcond,'S',sconds,'R134a');
    hconda = (hconds - ho)/eta + ho;
    hconda

    % work
    W = mdot*(hconda - ho);
    
    % work with no pressure loss
    so = py.CoolProp.CoolProp.PropsSI('S','P',Pi,'H',ho,'R134a');
    sconds = so;
    hconds = py.CoolProp.CoolProp.PropsSI('H','P',Pcond,'S',sconds,'R134a');
    hconda = (hconds - ho)/eta + ho;
    WP = mdot*(hconda - ho);

    % energy used
    dE = W*dt;

    % energy used no pressure loss
    dEP = WP*dt;
end






%% MAIN
% FITNESS FUNCTION - power and time
%f = 0.5t + 0.5w
Ttarget = -30+273;

% geometry conditions (from system design)
Vnos = 23.4/1000; %m^3
Lnos = 6*12*0.0254; %m
Dnos = (4*Vnos/pi/Lnos)^(1/2);

Lpipe = 1.5;
Hpipe = 1;
DpipeO = 0.5*0.0254;
Dpipei = DpipeO - 0.030*2*0.0254;
kpipe = 15;
D = Dnos - 1*0.0254*2;
n = sqrt(Lpipe^2-Hpipe^2)/(pi*D);
p = Hpipe/n;
inc = atan(p/(pi*D));

%corrected Vnos
Vnos = Vnos - pi/4*(DpipeO^2)*Lpipe;

Ltank = Lnos + 2*1*0.0254;
Dtanki = Dnos;
DtankO = Dnos + 0.25*0.0254*2;
ktank = 152;

Lins = 0.05;
Dins = DtankO + 2*0.0254*2;
kins = 0.033;

Rpp = 0.1;

% starting conditions
Tnos = 18+273; %C
mnos = 11; %kg
Tamb = 24+273; %C
mdot = 0.1;
Tlimit = -25+273;

% simulation conditions
dt = 2; %s, time step
dl = 0.01; %m, pipe flow integration length
tol = 0.75; %C, tempe`   rature convergence tolerance
tolq = 1; %kW, heat transfer rate convergence tolerance
cycleCount = 2;

% initial iteration guess
%Pevap = 126863.5-75000; %pa
Pevap = 126863.5-60000; %pa
Pcond = 226863.5+400000; %pa

% EVAPORATION
% calculate refrigeration states
[s,h,T] = RefrigerationCycle(Pevap,Pcond);
To = T(4);
Ti = T(3);
so = s(4);
si = s(3);
eta = 0.6;

% initialize time series and guesses
t = 0;
qnsin = 0.09;
%qnsout = 1200;
qnsout = 1000;
series = [0,0,0,0,0,0,0,0,0,0,0];
cycles = 0;
idleSeries = [];
activeSeries = [0,0,0,0,0,0,0,0,0,0,0];

% compute transient time series 
while Tnos > Ttarget
    % calculate heat transfer time step
    [Tnos,x,qin,qout,Pnos,dPfn,hon] = HeatTransfer(dt,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,DpipeO,Ti,To,si,so,Lpipe,Hpipe,Dpipei,kpipe,mdot,qnsin,qnsout,tol,tolq,dl,Pevap,D,n,p,inc);
    [W,WP,dE,dEP] = Power(Pcond,Pevap,dPfn,mdot,hon,eta,dt);

    % update counters and guesses
    qnsin = qin;
    qnsout = qout;
    t = t+dt;

    % log
    fprintf('t = %f, Tnos = %f C, x = %f, Pnos = %f psi\n', t, Tnos-273, x, Pnos*0.000145038);
    series = [series; t, Tnos-273, x, Pnos*0.000145038, W, WP, series(end,7)+dE, dE, series(end,9)+dEP, dEP, dPfn*0.000145038];
end
series(1,:) = [];

% compute cycling time series - 10 cycles
t = 0;
while cycles < cycleCount
    while Tnos < Tlimit
        % compute heat addition rate from ambient to NOS
        [Tnos,x,qin,Pnos] = IdleHeating(dt,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,DpipeO,Ti,To,si,so,Lpipe,Hpipe,Dpipei,kpipe,mdot,qnsin,qnsout,tol,tolq,dl,Pevap,D,n,p,inc);

        % update counters and guesses
        qnsin = qin;
        t = t+dt;
    
        % log
        fprintf('t = %f, Tnos = %f C, x = %f, Pnos = %f psi\n', t, Tnos-273, x, Pnos*0.000145038);
        idleSeries = [idleSeries; t, Tnos-273, x, Pnos*0.000145038,0,0,0,0,0,0,0];
    end
    while Tnos > Ttarget
        % calculate heat transfer time step
        [Tnos,x,qin,qout,Pnos,dPfn,hon] = HeatTransfer(dt,Tnos,Vnos,mnos,Lnos,Dtanki,Tamb,Lins,Dins,kins,Rpp,DtankO,ktank,DpipeO,Ti,To,si,so,Lpipe,Hpipe,Dpipei,kpipe,mdot,qnsin,qnsout,tol,tolq,dl,Pevap,D,n,p,inc);
        [W,WP,dE,dEP] = Power(Pcond,Pevap,dPfn,mdot,hon,eta,dt);
    
        % update counters and guesses
        qnsin = qin;
        qnsout = qout;
        t = t+dt;
    
        % log
        fprintf('t = %f, Tnos = %f C, x = %f, Pnos = %f psi\n', t, Tnos-273, x, Pnos*0.000145038);
        activeSeries = [activeSeries; t, Tnos-273, x, Pnos*0.000145038, W, WP, series(end,7)+dE, dE, series(end,9)+dEP, dEP, dPfn*0.000145038];
        activeSeries(1,:) = [];

    end
    cycleSeries = [idleSeries; activeSeries];
    idleSeries = [];
    activeSeries = [0,0,0,0,0,0,0,0,0,0,0];
    cycles = cycles + 1;
end


%% PLOTTING
% TRANSIENT
% Temperature time series
figure(1)
plot(series(:,1)/60,series(:,2),color='b')
title('NOS Temperature Time Series')
xlabel('Time (min)')
ylabel(['NOS Temperature (' char(176) 'C)']);

% Pressure time series
figure(2)
plot(series(:,1)/60,series(:,4),color='c')
title('NOS Pressure Time Series')
xlabel('Time (min)')
ylabel('NOS Pressure (psi)');

% Quality time series
figure(3)
plot(series(:,1)/60,series(:,3),color='k')
title('NOS Quality Time Series')
xlabel('Time (min)')
ylabel('NOS Quality');

% Work Time Series
figure(4)
hold on
plot(series(:,1)/60,smooth(series(:,5), 0.1, 'lowess'),color='g')
plot(series(:,1)/60,smooth(series(:,6), 0.1, 'lowess'),color='b')
title('Compressor Power Draw Time Series')
xlabel('Time (min)')
ylabel('Compressor Power Draw (kW)');
legend('With Pressure Loss', 'Without Pressure Loss')
hold off

% Energy Time Series
figure(5)
hold on
plot(series(:,1)/60,smooth(series(:,7), 0.1, 'lowess'),color='g')
plot(series(:,1)/60,smooth(series(:,9), 0.1, 'lowess'),color='b')
title('Total Energy Time Series')
xlabel('Time (min)')
ylabel('Total Energy (kJ)');
legend('With Pressure Loss', 'Without Pressure Loss')
hold off

% Pressure Loss Time Series
figure(6)
plot(series(:,1)/60,smooth(series(:,11), 0.1, 'lowess'),color='m')
title('Pressure Loss Time Series')
xlabel('Time (min)')
ylabel('Pressure Loss (psi)');

% CYLCLING
% Temperature time series
figure(1)
plot(cycleSeries(:,1)/60,cycleSeries(:,2),color='b')
title('NOS Temperature Time Series')
xlabel('Time (min)')
ylabel(['NOS Temperature (' char(176) 'C)']);

% Pressure time series
figure(2)
plot(cycleSeries(:,1)/60,cycleSeries(:,4),color='c')
title('NOS Pressure Time Series')
xlabel('Time (min)')
ylabel('NOS Pressure (psi)');

% Quality time series
figure(3)
plot(cycleSeries(:,1)/60,cycleSeries(:,3),color='k')
title('NOS Quality Time Series')
xlabel('Time (min)')
ylabel('NOS Quality');

% Work Time Series
figure(4)
hold on
plot(cycleSeries(:,1)/60,smooth(cycleSeries(:,5), 0.1, 'lowess'),color='g')
plot(cycleSeries(:,1)/60,smooth(cycleSeries(:,6), 0.1, 'lowess'),color='b')
title('Compressor Power Draw Time Series')
xlabel('Time (min)')
ylabel('Compressor Power Draw (kW)');
legend('With Pressure Loss', 'Without Pressure Loss')
hold off

% Energy Time Series
figure(5)
hold on
plot(cycleSeries(:,1)/60,smooth(cycleSeries(:,7), 0.1, 'lowess'),color='g')
plot(cycleSeries(:,1)/60,smooth(cycleSeries(:,9), 0.1, 'lowess'),color='b')
title('Total Energy Time Series')
xlabel('Time (min)')
ylabel('Total Energy (kJ)');
legend('With Pressure Loss', 'Without Pressure Loss')
hold off

% Pressure Loss Time Series
figure(6)
plot(cycleSeries(:,1)/60,smooth(cycleSeries(:,11), 0.1, 'lowess'),color='m')
title('Pressure Loss Time Series')
xlabel('Time (min)')
ylabel('Pressure Loss (psi)');

%% TODO
% check refrigerant convection - mid values (need integral?)
    % did lengthwise integral with constant surface temp
    % solved for surface temp for a guess q
    % see HeatExchanger_Flow file, comments in Main for implementation
    % questionable, and run time is very bad
    % yes you need to integrate on length, look at boiling models in text
% update refrigerant convection with pressure drop - flow model, then int
    % nah, assume constant that would suck
% update vertical cylinder assumption if needed - assumed thin BL
% update for correct refrigerant
% update refrigerant flow models file - check mdot selection
% update for nitrous data - REFPROP through U of C?
    % maybe, ask Aggrey
% optimize runtime - coolprop calls no bueno, verrrrry slow
% correct refrigeration cycle for pressure drop - recalc states
% add transient zone for reynolds number - currently Re<2300 or Re>10000
% advice on variable mass flow rate - changing ref. cycle each time step
% correct displacement volume of nitrous due to the evaporator - at start

% look at boiling condensing HT models in textbook 