% the main.m
clear all
algorithmDir_basic = 'DE_prox';
algorithmDir = 'DE_prox';
optimalChart_interval = 100;

%%
% Global parameters setting

%CEC2005 include 25 functions
%CEC2011 include 22 functions
%CEC2013 include 28 functions
%CEC2014 include 30 functions
%CEC2017 include 30 functions
benchmark = 2017;
problemIndex = (1:30);
runTime = 51;
D = 30;
popuSize = 100;
FES = 1 * 10^4 * D;
addpath('../../');
%% For CEC2011
if benchmark ~= 2011
    pathOfBoxPlot = ['./',algorithmDir,'_CEC',num2str(benchmark),'_D',num2str(D),'_Mean&Std_and_Box-Plot.xls'];
    pathOfConvergence = ['./',algorithmDir,'_CEC',num2str(benchmark),'_D',num2str(D),'_Convergence.xls'];
else
    %% init table for CEC2011
    pathOfBoxPlot = ['./',algorithmDir,'_CEC',num2str(benchmark),'_Mean&Std_and_Box-Plot.xls'];
    pathOfConvergence = ['./',algorithmDir,'_CEC',num2str(benchmark),'_Convergence.xls'];
end
%% init table
tableOfBoxPlot = initTableforBoxPlot(length(problemIndex),runTime);
tableOfConvergence = initTableforConvergence(length(problemIndex));
rand('seed', sum(100 * clock));
%rng('shuffle');

val_2_reach = 10^(-8);
for func = 30%problemIndex
    resultChart = [];
    f = func;
    fnum = f;
    optimum = func * 100.0;
    if benchmark == 2011
        if f == 11; fnum = 111; end
        if f == 12; fnum = 112; end
        if f == 13; fnum = 113; end
        if f == 14; fnum = 114; end
        if f == 15; fnum = 115; end
        if f == 16; fnum = 116; end
        if f == 17; fnum = 117; end
        if f == 18; fnum = 118; end
        if f == 19; fnum = 119; end
        if f == 20; fnum = 1110; end
        if f == 21; fnum = 12; end
        if f == 22; fnum = 13; end
        CEC2011Dir = 'CEC/CEC_2011/CEC 2011 Benchmarks';
        addpath(['../../', CEC2011Dir]);
        GetFunInfoAddr = @CEC2011get_fun_info;
        [LB, UB, dem, opt_f, err] = GetFunInfoAddr(fnum);
        D = dem;disp(D);
        %% ============================================================================================================================================================
        %% ============================================================================================================================================================
        popuSize = 100;  %fix this line when running CEC_2011
        FES = 1 * 10^4 * D; %fix this line when running CEC_2011
        %% ============================================================================================================================================================
        %% ============================================================================================================================================================
    end
    for run_id = 1 : runTime
        %% init of each run
        % FES count
        nFES = 0;
        nfes = nFES;
        % iteration count
        iter = 0;
        optimalChart = [];
        % Population initialization
        [popold, lu, rgo, o, A, M, a, alpha, b] = Initialization(popuSize, D, benchmark, f);
        % dimension for certain problem
        D = size(popold,2);
        % program start here
        if rgo ~= 0
            rgo = rgo(f);
        end
        bsf_fit_var = 1e+30;
        % Population evaluation
        [popuFitness] = Evaluation(popold, benchmark, f);
        
        fit = popuFitness;
        fit = fit';
        F = 0.5;
        CR = 0.9;
        popu= popold;
        
        % trajectory map
        if D == 2
            mkdir traj
            imgPath = ['traj_' num2str(benchmark) '_' num2str(problemIndex)];
            mkdir(imgPath);
        end
        
        % get bsf_fit_var
        for i = 1 : popuSize
            nFES = nFES + 1;
            if fit(i) < bsf_fit_var
                bsf_fit_var = fit(i);
            end
            if mod(nFES, optimalChart_interval) == 0
                optimalChart = [optimalChart;bsf_fit_var];
            else
                if nFES == FES
                    optimalChart = [optimalChart;bsf_fit_var];
                end
            end
            if nFES > FES; break; end
        end
        
        g_num = 0;
        % Main loop
        while nFES < FES
            
            p = nFES/FES;
            p = 0.5*(1-p);
            
            
                % 获取前一半
                [AllFitnessSorted,IndexSorted] = sort(fit);
                N = popuSize;
                N_half = round(N/2);
                index = floor((1-p)*N_half);
                
                goodPops=popu(IndexSorted(1:N_half),:);
                worstPops=popu(IndexSorted(N_half+1:N),:);
                
                a = randperm(N_half);
                b = randperm(N_half);
                
                gP = goodPops(a,:);
                gp1 = gP(1:index,:);
                gp2 = gP(index+1:N_half,:);
                
                wP = worstPops(b,:);
                wp1 = wP(1:N_half-index,:);
                wp2 = wP(N_half-index+1:N_half,:);
                
                V_pops = [gp1;wp1];

                U = [gp2;wp2];
                    % Get indices for mutation
                [r1, r2, r3] = getindex(N_half);
                     % Implement DE/rand/1 mutation
                V = V_pops(r1, :) + F * (V_pops(r2, :) - V_pops(r3, :));
                    % Check whether the mutant vector violates the boundaries or not
                [V] = BoundaryDetection(V,lu);
                    % Implement binomial crossover
                for i = 1:N_half
                    j_rand = floor(rand * D) + 1;
                    t = rand(1, D) < CR;
                    t(1, j_rand) = 1;
                    t_ = 1 - t;
                    U(i, :) = t .* V(i, :) + t_ .* U(i, :);
                    %X(IndexSorted(i), :) = t .* V(i, :) + t_ .* X(IndexSorted(i), :);
                end

                % evaluate population
                [fit_U] = Evaluation(U, benchmark, f); % 较差个体的fitness
                fit_U=fit_U'; 

                % get bsf_fit_var
                for i = 1:N_half
                    nFES = nFES + 1;
                    if fit_U(i) < bsf_fit_var
                        bsf_fit_var = fit_U(i);
                    end
                    if mod(nFES, optimalChart_interval) == 0
                        optimalChart = [optimalChart;bsf_fit_var];
                    else
                        if nFES == FES
                            optimalChart = [optimalChart;bsf_fit_var];
                        end
                    end
                    popu_num = IndexSorted(N_half+i); % 差种群个体的编号
                    if fit_U(i, :) <= fit(popu_num, :)
                        popu(popu_num, :) = U(i, :);
                        fit(popu_num, :) = fit_U(i, :);
                    end

                    if nFES > FES; break; end
                end
                
            
            %% Survival
            % For minimization problems
            optimal = bsf_fit_var;
            fprintf('%s problem %5.0f time %5.0f |%5.0f %5.5f-----> %9.16f\n', algorithmDir, f,run_id,nFES,index,optimal);
            if D == 2
                showTraj(f, nfes, popu(1:popuSize,:), optimal, optimal, imgPath);
            end
        end
%         bsf_error_val = optimalChart - rgo;
%         if bsf_error_val < val_2_reach
%             bsf_error_val = 0;
%         end
%         optimalChart = bsf_error_val;
        resultChart = [resultChart,optimalChart];
    end
    boxPlotChart = resultChart(end,:);
    convergenceChart = mean(resultChart,2);
    if benchmark ~= 2011
        path = ['./',algorithmDir_basic,'_',algorithmDir,'_CEC',num2str(benchmark),'_D',num2str(D),'.xls'];
    else
        path = ['./',algorithmDir_basic,'_',algorithmDir,'_CEC',num2str(benchmark),'.xls'];
    end
    sheetName = ['CEC_' num2str(benchmark) '_F'];
    xlswrite(path,resultChart,[sheetName,num2str(f),'_D',num2str(D)]);
    
    %% save result
    tableOfBoxPlot = updateTableforBoxPlot(tableOfBoxPlot,boxPlotChart,pathOfBoxPlot,rgo,f);
    tableOfConvergence = updateTableforConvergence(tableOfConvergence,convergenceChart,pathOfConvergence,rgo,f);
end

function tableB = initTableforBoxPlot(problemNum, runNum)
row = 3 + problemNum;
col = 1 + runNum;
tableB = cell(row, col);
tableB{1,2} = 'mean';
tableB{1,3} = 'std';
for i = 1:problemNum
    tableB{i+1,1} = ['F' num2str(i)];
end
end

function tableC = initTableforConvergence(problemNum)
col = problemNum;
tableC = cell(1,col);
for i = 1:col
    tableC{1,i} = ['F' num2str(i)];
end
end

function tableB = updateTableforBoxPlot(tableB,data,path,rgo,problemIndex)
[row,col] = size(data);
% rgo = repmat(rgo(problemIndex),1,col);
rgo = repmat(rgo,1,col);
data = data - rgo;
dataMean = mean(data,2);
dataStd = std(data,0,2);
Dcol = ones(1,col);
Drow = ones(1,row);
dataMean = mat2cell(dataMean,Drow,[1]);
dataStd = mat2cell(dataStd,Drow,[1]);
tableB(problemIndex+1,2) = dataMean;
tableB(problemIndex+1,3) = dataStd;
dataCol = mat2cell(data,[1],Dcol);
tableB(problemIndex+1,4:3+col) = dataCol;
xlswrite(path,tableB);
end

function tableC = updateTableforConvergence(tableC,data,path,rgo,problemIndex)
[row,col] = size(data);
% rgo = repmat(rgo(problemIndex)',row,1);
rgo = repmat(rgo',1,col);
data = data - rgo;
Drow = ones(1,row);
dataCol = mat2cell(data,Drow,[1]);
tableC(2:row+1,problemIndex) = dataCol;
xlswrite(path,tableC);
end

function showTraj(problemIndex, nFES, popu, pBest, optimal, imgPath)
popuSize = size(popu,1);
func_plot(problemIndex);
h = plot(popu(:,1),popu(:,2),'.r','MarkerSize',15);
% r = sectionStep;
% for i = 1:popuSize
% x = popu(i,1);
% y = popu(i,2);
% rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1),axis equal
hold on
% plot(optimal(1,1),optimal(1,2),'-rp','MarkerSize',20);
% plot(pBest(:,1),pBest(:,2),'.b','MarkerSize',15);
if mod(nFES / popuSize, 1) == 0
    saveas(gcf,[[imgPath '/'],num2str(nFES / popuSize) '.png']);
    saveas(gcf,[[imgPath '/'],num2str(nFES / popuSize) '.eps'], 'epsc');
end
hold off
end
