% Brock Lumbers Thesis 24.08.2020
% Mathematical Modelling and Simulation of Catalyst Deactivation For The Continuous Thermo-catalytic Decomposition of Methane

% Calculation of methane equilibrium conversion

% Please note that execution times may be up to 15 seconds

R = 8.31447;    % J/mol K
P = 1;          % Bar
P2 = 3;         % Bar
P3 = 6;         % Bar
P4 = 9;         % Bar
T0 = 603;       % K
Tmax = 1403;    % K
stpsze = 10;    % stepsize
Xeq = 0;	    % initial conversion
Xeq2 = 0;	    % initial conversion
Xeq3 = 0;	    % initial conversion
Xeq4 = 0;	    % initial conversion
 
Trange = [(T0-stpsze):stpsze:Tmax] ;   % Temperature Range
 
for i = 0:stpsze:(Tmax - T0)
    
        Xprev = Xeq;
        Xprev2 = Xeq2;
        Xprev3 = Xeq3;
        Xprev4 = Xeq4;
        
        T = T0+i ; 

        % Villacampa 2003 Change in Gibbs Energy Correlation
        G = 89658.88 - 102.27*(T) - 0.00428*(T)^(2) - 2499358.99/(T) ;

        % Calculation of Equilibrium Constant
        K =  (exp((-1*G)/(R*T))); 
        
        % Solving for the Equilibrium Methane Conversion
        syms x x2 x3 x4 
        assume(x >= 0)
        assume(x2 >= 0)
        assume(x3 >= 0)
        assume(x4 >= 0)
        Xe = solve(((P)*((2*x)/((1+2*x))^2))/(((1-x)/(1+2*x)))==K,x) ;
        Xe2 = solve(((P2)*((2*x2)/((1+2*x2))^2))/(((1-x2)/(1+2*x2)))==K,x2) ;
        Xe3 = solve(((P3)*((2*x3)/((1+2*x3))^2))/(((1-x3)/(1+2*x3)))==K,x3) ;
        Xe4 = solve(((P4)*((2*x4)/((1+2*x4))^2))/(((1-x4)/(1+2*x4)))==K,x4) ;
        Xeq = [Xprev] ;
        Xeq2 = [Xprev2] ;
        Xeq3 = [Xprev3] ;
        Xeq4 = [Xprev4] ;
       
       Xeq = [Xeq Xe] ;
       Xeq2 = [Xeq2 Xe2] ;
       Xeq3 = [Xeq3 Xe3] ;
       Xeq4 = [Xeq4 Xe4] ;
end
 
plot(Trange,Xeq,'b');
hold on
plot(Trange,Xeq2,'b');
hold on
plot(Trange,Xeq3,'b');
hold on
plot(Trange,Xeq4,'b');
axis([T0,Tmax,0,1]);
xlabel('Temperature (K)');
ylabel('Methane Conversion (-)');
text(Trange(27), Xeq(27), '1 bar', 'BackgroundColor', 'w', 'rotation', 72);
text(Trange(33), Xeq2(33), '3 bar', 'BackgroundColor', 'w', 'rotation', 69);
text(Trange(37), Xeq3(37), '6 bar', 'BackgroundColor', 'w', 'rotation', 68);
text(Trange(39), Xeq4(39), '9 bar', 'BackgroundColor', 'w', 'rotation', 67);
