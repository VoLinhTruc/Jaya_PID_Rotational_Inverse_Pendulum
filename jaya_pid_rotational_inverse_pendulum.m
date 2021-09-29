clc;
clear;


Beq = 0.0035842;
g = 9.81;
Jeq = 0.004;
Jm = 67.7*10^-7;
Kg = 102/16;
Km = 36.4*10^-3;
L = 10*10^-2;
m = 50*10^-3;
r = 10*10^-2;
Rm = 1.11;
ng = 0.9;

a = Jeq + m*r^2 + Jm*Kg^2*ng;
b = m*L*r;
c = (4/3)*m*L^2;
d = m*L*g;
e = Beq + (ng*Kg^2*Km^2)/Rm;
f = Km*Kg*ng/Rm;

A = [0 0 1 0; 0 0 0 1; 0 (b*d)/(a*c-b^2) -(c*e)/(a*c-b^2) 0; 0 (a*d)/(a*c-b^2) -(b*e)/(a*c-b^2) 0];
B = [0 0 c*f/(a*c-b^2) b*f/(a*c-b^2)]';
C = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
D = [0; 0; 0; 0];



global P_alpha;
global I_alpha;
global D_alpha;

global P_theta;
global I_theta;
global D_theta;

Jaya();

function Jaya()

    pop = 10; % population size
    var = 6; % no. of design variables
    maxGen = 50;
    mini = 0 * ones(1, var);
    maxi = 50 * ones(1, var);
    jaya_table = zeros(pop, var);
    
    RUNS = 1;
    runs = 0;

    while (runs < RUNS)

        % initialize jaya_table randomly
        for i = 1:var % P_alpha and P_theta
            if (i == 1)
                jaya_table(:, i) = 25 + (mini(i) + (maxi(i) - mini(i)) * rand(pop, 1))/50;
            elseif((i == 2) || (i == 5))
                jaya_table(:, i) = (mini(i) + (maxi(i) - mini(i)) * rand(pop, 1))/100;
            else
                jaya_table(:, i) = (mini(i) + (maxi(i) - mini(i)) * rand(pop, 1))/50;                
            end
        end

        gen = 0;
        f = pidSimulationRun(jaya_table);

        while (gen < maxGen)
            jaya_table_new = updatePopulation(jaya_table, f);
            jaya_table_new = trimr(mini, maxi, jaya_table_new);
            f_new = pidSimulationRun(jaya_table_new);

            % after updating, keep the optimizer f value
            for i = 1:pop
                if (f_new(i) < f(i))
                    jaya_table(i, :) = jaya_table_new(i, :);
                    f(i) = f_new(i);
                end
            end

            %disp('%%%%%%%% Final population %%%%%%%%%');
            %disp([jaya_table, f]);

            % if norm(jaya_table - jaya_table_new) < 0.001
            %     break;
            % end

            f_new = []; jaya_table_new = [];
            gen = gen + 1;
            f_optimized(gen) = min(f);
        
            %{
            clf;
            hold on;
            plot(0, 0, '.');
            plot(pop + 1, 0, '.');
            for i = 1:pop
                plot(i, jaya_table(i, 1), '.');
                ylim([-5 5]);
                pause(0.05);
            end
            hold off;
            
            disp(gen)
            %}
            
        end

        runs = runs + 1;
        [val, ind] = min(f_optimized);
        %Fes(runs) = pop * ind;
        best(runs) = val;
        [f_min, f_min_idx] = min(f);
        x_opt(runs, :) = jaya_table(f_min_idx, :);
    end

    bbest = min(best);
    mbest = mean(best);
    wbest = max(best);
    %stdbest = std(best);
    %mFes = mean(Fes);
    %stdFes = std(Fes);

    fprintf('\n runs=%f', runs);
    fprintf('\n best=%f', bbest);
    fprintf('\n mean=%f', mbest);
    fprintf('\n worst=%f', wbest);
    %fprintf('\n std. dev.=%f', stdbest);
    %fprintf('\n mean function evaluations=%f', mFes);
    disp(' ');
    disp('x_opt: ');
    disp(x_opt);
end

function [z] = trimr(mini, maxi, jaya_table)
    [row, col] = size(jaya_table);

    for i = 1:col
        jaya_table(jaya_table(:, i) < mini(i), i) = mini(i);
        jaya_table(jaya_table(:, i) > maxi(i), i) = maxi(i);
    end

    z = jaya_table;
end

function [jaya_table_new] = updatePopulation(jaya_table, f)
    [row, col] = size(jaya_table);
    % ---------------------------------
    [b, bIndex] = min(f);
    [w, wIndex] = max(f);
    % ---------------------------------
    best = jaya_table(bIndex, :);
    worst = jaya_table(wIndex, :);
    jaya_table_new = zeros(row, col);

    for i = 1:row
        for j = 1:col
            r = rand(1, 2);
            jaya_table_new(i, j) = jaya_table(i, j) + r(1)*(best(j) - abs(jaya_table(i, j))) - r(2)*(worst(j) - abs(jaya_table(i, j)));
        end
    end

end

function [f] = pidSimulationRun(jaya_table)
    global P_alpha;
    global I_alpha;
    global D_alpha;

    global P_theta;
    global I_theta;
    global D_theta;
    
    [r, c] = size(jaya_table);
    for i = 1:r
        P_alpha = jaya_table(i, 1);
        I_alpha = jaya_table(i, 2);
        D_alpha = jaya_table(i, 3);
        P_theta = jaya_table(i, 4);
        I_theta = jaya_table(i, 5);
        D_theta = jaya_table(i, 6);
        
        sim('pid_inverse_pendulum');
        
        z(i) = norm_evaluation.data;
    end
    
    f =  z';
    %{
    fprintf('\n max_x_dot =%f', max(abs(x_dot.data)));
    fprintf('\n P=%f \t I=%f \t D=%f \n', P, I, D);
    disp(f)
    %}
end

function [f] = myobj(jaya_table)
    [r, c] = size(jaya_table);

    for i = 1:r
        for j = 1:c - 2
            y = jaya_table(i, j)^2 + jaya_table(i, j + 1)^2 + jaya_table(i, j + 2)^2 + 5;
        end
        z(i) = y;
    end
    
    f = z';
end
