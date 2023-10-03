function metabolism_model()
 
    Vb = 50;  % Volume of distribution for parent drug
    Vm = 50;  % Volume of distribution for metabolite 
    k = 10;   % Rate constant for elimination of parent drug
    km = 20; % Rate constant for elimination of metabolite
    fc = 0.3; % Fraction of metabolite converted back to parent drug
    fm = 0.2; % Fraction of parent drug converted to metabolite

    % Simulation parameters
    tspan = [0 60];  % Simulation time (minutes)
    initial_conditions = [1000 0];  % Initial drug mass from IV bolus injection
    
    % Solve the differential equations using ode45
    [t, drug_amounts] = ode45(@(t, y) ode_equations(y, Vb, Vm, k, km, fc, fm), tspan, initial_conditions);

    % Plot drug amounts in the central and peripheral compartments
    figure;
    plot(t, drug_amounts(:, 1), 'r-', t, drug_amounts(:, 2), 'b-');
    xlabel('Time');
    ylabel('Concentration');
    legend('Parent Drug', 'Metabolite');
    title('Concentration vs. Time (Metabolism and Excretion)');

end

function dydt = ode_equations(y, Vb, Vm, k, km, fc, fm)
   
    C = y(1) / Vb;  
    Cm = y(2) / Vm; 

    dydt = [
             -k*C + fc*km*Cm*(Vm/Vb);
             -km*Cm + fm*k*C*(Vb/Vm);
            ]; 
end