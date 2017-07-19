%% CH E 454 - Lab 2D Data File
clc
clear
close all

%% Load Raw Data

load('CHE454_2DRaw')
Run{1} = M;
Run{2} = T;
Run{3} = W;
Run{4} = F;
Run{5} = S;

%% Constants

R = 8.314; %J/mol/K

%Experimental
V = 0.637; %L
dt = 5; %seconds
C_AT = 0.05; %mol/L
C_BT = 0.0576; %mol/L

%Literature
Ea = [54448 47100.00 38867.35 45600.00 47861.00];
A  = [1.82E8 2.9000E7 7.1255E5 3.7400E7 2.7466E7];

%Steady-State Time

startSS = [
    (259-18+1),0,(925-684+1); %M
    (219-18+1),0,(1148-838+1); %T
    (419-18+1),0,(1338-1137+1); %W
    (319-18+1),0,(1290-989+1); %F
    (260-18+1),0,0]; %S

endSS = [
    (347-18+1),0,(999-684+1); %M
    (412-18+1),0,(1309-838+1); %T
    (603-18+1),0,(1631-1137+1); %W
    (572-18+1),0,(1519-989+1); %F
    (534-18+1),0,0,]; %S



%% Data Structure Development

%pre-allocation

time = cell(length(Run),length(Run{1}));
temp = cell(length(Run),length(Run{1}));
flowA = cell(length(Run),length(Run{1}));
flowB = cell(length(Run),length(Run{1}));
cond = cell(length(Run),length(Run{1}));

limit = 1E-6*[0.2 1 1 0.1 1];
%cellArray{run,part}
%for parts: 1 = flow, 2 = RTD, 3 = temp

for i = 1:length(Run)
    for j = 1:length(M)
        if ~isempty(Run{1,i}{1,j})
            time{i,j} = Run{1,i}{1,j}(:,1);
            temp{i,j} = Run{1,i}{1,j}(:,2);
            flowA{i,j} = Run{1,i}{1,j}(:,3);
            flowB{i,j} = Run{1,i}{1,j}(:,4);
            cond{i,j} = Run{1,i}{1,j}(:,5);
        else
            time{i,j} = [];
            temp{i,j} = [];
            flowA{i,j} = [];
            flowB{i,j} = [];
            cond{i,j} = [];
        end
    end
end

%% Calculations

for i = 5:-1:1
    
    for j = [2,1,3]
        %% Preliminary
        
        if i == 5 && j == 3
            break
        end
        
        %temporary variables
        xtime = time{i,j};
        xtemp = temp{i,j};
        xflowA = flowA{i,j};
        xflowB = flowB{i,j};
        xcond = cond{i,j};
        
        %Steady State
        
        %flow and temp
        totalFlow{i,j} = xflowA + xflowB;
        xtotalFlow = totalFlow{i,j};
        
        temp_av{i,j} = mean(temp{i,j});
        xtemp_av = temp_av{i,j};
        
        totalFlow_av{i,j} = mean(xtotalFlow);
        xtotalFlow_av = totalFlow_av{i,j};
        
        %space time
        tao{i,j} = 1000*V*60 ./ mean(xflowA);
        xtao = tao{i,j};
        
        %% Residence Time Distribution
        if j == 2
            
            %concentration
            C_A{i,j} = xcond ./(192.5.*(1+0.01722.*(xtemp-20)));
            intC_A{i} = homer(C_A{i,j},dt);
            
            %RTD
            RTD_exp{i} = C_A{i,j}/intC_A{i};
            RTD_pm{i,j} = exp(-xtime./xtao)./xtao;
            
            %curve fitting

            
            RTD_fit{i,j} = fit(xtime,RTD_exp{i},'linearint');
            RTD_act{i,j} = RTD_fit{i,j}(xtime);
            dRTD_act{i,j} = gradient(RTD_act{i,j});
            cutStart = find(max(RTD_act{i,j})==RTD_act{i,j});
            cutEnd = length(xtime);
            shift = min(RTD_act{i,j});
            RTD_cut{i,j} = RTD_act{i,j}(cutStart:cutEnd);
            adRTD_cut{i,j} = abs(gradient(RTD_cut{i,j}));
            adRTD_cut_fit{i,j} = fit(xtime(cutStart:cutEnd),adRTD_cut{i,j},...
                'exp1');
            figure
            plot(adRTD_cut_fit{i,j},xtime(cutStart:cutEnd),adRTD_cut{i,j})
            title(i)
            adRTD_cut{i,j} = adRTD_cut_fit{i,j}(xtime(cutStart:cutEnd));
            cutEnd = min(find(adRTD_cut{i,j}<limit(i)));
            if isempty(cutEnd)
                cutEnd = length(xtime);
            end
            RTD_cut{i,j} = RTD_act{i,j}(cutStart:cutEnd);
            RTD_shift{i,j} = RTD_act{i,j} - shift;
            RTD_shiftcut{i,j} = RTD_cut{i,j} - shift;
            
            dif = length(RTD_cut{i,j}) - length(RTD_act{i,j});
            pad = zeros(-dif,1);
            xRTD_cut_pad = [RTD_cut{i,j}; pad];
            
            xRTD_shiftcut_pad = [RTD_cut{i,j}; pad];
            EXPORT_RTD{i} = [xtime RTD_pm{i,j}, RTD_act{i,j}, xRTD_cut_pad, RTD_shift{i,j},...
                xRTD_shiftcut_pad];
            
                    
                
            
            %mean residence time
            
            xcutTime = xtime(1:(cutEnd-cutStart+1));
            tm(i) = homer(xtime.*RTD_act{i,j},dt);
            tm_cut(i) = homer(xcutTime.*RTD_cut{i,j},dt);
            tm_shift(i) = homer(xtime.*RTD_shift{i,j},dt);
            tm_shiftcut(i) = homer(xcutTime.*RTD_shiftcut{i,j},dt);
            
            
            %plots
            figure
            plot(RTD_fit{i,j}, xtime, RTD_pm{i,j})
            title(i)
            
            
            %export
            taoArray(i,j) = xtao;
            
            
        else
            
            %% Experimental Conversion
            
            
            
            %concentration
            C_A0{i,j} = xflowB * C_AT ./ xtotalFlow;
            xC_A0 = C_A0{i,j};
            
            C_B0{i,j} = xflowA * C_BT ./ xtotalFlow;
            xC_B0 = C_B0{i,j};
            
            expC_Ct = ((192.5 * ( 1 + 0.01722 * (xtemp - 20)) .* xC_A0)...
                - xcond) ./ (123.82 + 1.6651564 .* (xtemp-20));
            
            xC_A = xC_A0 - expC_Ct;
            xC_B = xC_B0 - expC_Ct;
            
            %conversion
            ExpConv{i,j} = 1 - xC_A ./ xC_A0;
            xExpConv = ExpConv{i,j};
            
            %k
            k_exp{i,j} = xtotalFlow.*(xC_A0 - xC_A) ./ ...
                (V.*1000.*xC_A.*xC_B*60);
            xk_exp = k_exp{i,j};
            
            %% Literature
            
            %Steady State
            
            start = startSS(i,j);
            stop = endSS(i,j);
            
            if cutEnd - cutStart < endSS(i,j) - startSS(i,j)
                cutEnd = cutStart + (stop - start);
                RTD_cut{i,2} = RTD_act{i,2}(cutStart:cutEnd);
                RTD_shift{i,2} = RTD_act{i,2} - shift;
                RTD_shiftcut{i,2} = RTD_cut{i,2} - shift;
            end
            
            %time
            sstime{i,j} = 5*0:5:5*(stop-start);
            xsstime = sstime{i,j};
            
            %conductivity
            sscond{i,j} = xcond(start:stop);
            
            %experimental conversion
            ssExpConv{i,j} = xExpConv(start:stop);
            
            %inlet concentration
            ssC_A0{i,j} = xC_A0(start:stop);
            ssC_B0{i,j} = xC_B0(start:stop);
            
            xssC_A0 = ssC_A0{i,j};
            xssC_B0 = ssC_B0{i,j};
            
            %RTD - exp and pm (RTD adjusted based on cutoff point and time)
            
            %pm
            RTD_pm{i,j} = exp(-xsstime./xtao)./xtao;
            
            
            %experimental
            
            cutoff = length(xsstime);
            %for part 1 of run i
            if j == 1
                RTD_fit{i,j} = RTD_act{i,2}(1:cutoff);
                RTD_shift{i,j} = RTD_shift{i,2}(1:cutoff);
                RTD_cut{i,j} = RTD_cut{i,2}(1:cutoff);
                RTD_shiftcut{i,j} = RTD_shiftcut{i,2}(1:cutoff);
            end
            
            if j == 3
                RTD_fit{i,j} = RTD_act{i,2}(1:cutoff);
                RTD_shift{i,j} = RTD_shift{i,2}(1:cutoff);
                RTD_cut{i,j} = RTD_cut{i,2}(1:cutoff);
                RTD_shiftcut{i,j} = RTD_shiftcut{i,2}(1:cutoff);
            end
            
            for it = 1:length(xsstime)
                F_fit{i,j}(it) = homer(RTD_fit{i,j}(1:it),dt);
                F_fit_shift{i,j}(it) = homer(RTD_shift{i,j}(1:it),dt);
                F_pm{i,j}(it) = homer(RTD_pm{i,j}(1:it),dt);
            end
            
            for it = 1:length(xsstime)
                F_fit_shiftcut{i,j}(it) = homer(RTD_shiftcut{i,j}(1:it),dt);
                F_fit_cut{i,j}(it) = homer(RTD_cut{i,j}(1:it),dt);
            end
            
            for s = 1:length(Ea)
                
                k_lit{i,j,s} = A(1,s).*exp(-Ea(1,s)./R./(xtemp_av + 273));
                xk_lit = k_lit{i,j,s};
                
                %Perfectly Mixed Model
                for it = 1:length(xC_A0)
                    a = 1;
                    b = -(xC_A0(it)+xC_B0(it)+1/xk_lit/xtao);
                    pmB{i,j,s}(it,1) = b;
                    c = xC_A0(it)*xC_B0(it);
                    pmC_C = (-b - sqrt(b^2-4*a*c))/2/a;
                    pmC_Ct{i,j,s}(it,1) = pmC_C;
                    pmC_A{i,j,s}(it,1) = xC_A0(it) ...
                        - pmC_C;%second (-) root chosen
                    Conv_pm{i,j,s}(1,it) = 1 -...
                        pmC_A{i,j,s}(it,1)./xC_A0(it);
                end
                ssConv_pm{i,j,s} = Conv_pm{i,j,s}(start:stop);
                
                %Segregated Model
                sA = exp(xk_lit.*xsstime.*(xssC_B0'-xssC_A0'));
                Conv_s{i,j,s} = (1-sA) ./ (xssC_A0'./xssC_B0' - sA);
                
                if length(xsstime) > cutStart - cutEnd
                    Conv_s_cut{i,j,s} = Conv_s{i,j,s};
                else
                    Conv_s_cut{i,j,s} = Conv_s{i,j,s}(cutStart:cutEnd);
                end
                
                %Max-Mix Model
                Conv_mm{i,j,s}(1)=0;
                
                for it = 2:length(sstime{i,j})
                    mmA = RTD_fit{i,j}(it)/(1-F_fit{i,j}(it));
                    mmA_shift = RTD_shift{i,j}(it)...
                        /(1-F_fit_shift{i,j}(it));
                    mmA_pm = RTD_pm{i,j}(it)...
                        /(1-F_pm{i,j}(it));
                    mmB = xk_lit*xssC_A0(it);
                    mmX = Conv_mm{i,j,s}(it-1);
                    Conv_mm{i,j,s}(it,1) = mmX - ...
                        dt*(mmA*mmX - mmB*(1-mmX)^2);
                    Conv_mm_shift{i,j,s}(it,1) = mmX - ...
                        dt*(mmA_shift*mmX - mmB*(1-mmX)^2);
                    Conv_mm_pm{i,j,s}(it,1) = mmX - ...
                        dt*(mmA_pm*mmX - mmB*(1-mmX)^2);

                end
                
                for it = 2:length(xsstime)
                    mmA_cut = RTD_cut{i,j}(it)/(1-F_fit_cut{i,j}(it));
                    mmA_shiftcut = RTD_shiftcut{i,j}(it)...
                        /(1-F_fit_shiftcut{i,j}(it));
                    Conv_mm_cut{i,j,s}(it,1) = mmX - ...
                        dt*(mmA_cut*mmX - mmB*(1-mmX)^2);
                    Conv_mm_shiftcut{i,j,s}(it,1) = mmX - ...
                        dt*(mmA_shiftcut*mmX - mmB*(1-mmX)^2);  
                end
                
                %Average Conversion
                
                avConv_exp{s}(i,j) = homer(ssExpConv{i,j}.*RTD_fit{i,j},dt);%/...
                    %homer(RTD_fit{i,j}',dt);
                avConv_pm{s}(i,j) = homer(ssConv_pm{i,j,s}.*RTD_fit{i,j}',dt);%/...
                    %homer(RTD_fit{i,j}',dt);
                
                %actual - unadjusted
                avConv_s{s}(i,j) = homer(Conv_s{i,j,s}.*RTD_fit{i,j}',dt);%/...
                    %homer(RTD_fit{i,j}',dt);
                avConv_mm{s}(i,j) = homer(Conv_mm{i,j,s}.*RTD_fit{i,j},dt);%/...
                    %homer(RTD_fit{i,j}',dt);
                
                %perfectly mixed RTD for all
                avConv_exp_pm{s}(i,j) = homer(ssExpConv{i,j}.*RTD_pm{i,j}',dt);%/...
                    %homer(RTD_pm{i,j}',dt);
                avConv_s_pm{s}(i,j) = homer(Conv_s{i,j,s}.*RTD_pm{i,j},dt);%/...
                    %homer(RTD_pm{i,j}',dt);
                avConv_mm_pm{s}(i,j) = homer(Conv_mm_pm{i,j,s}.*RTD_pm{i,j}',dt);%/...
                    %homer(RTD_pm{i,j}',dt);
               
                
                %cut RTD
                avConv_exp_cut{s}(i,j) = homer(ssExpConv{i,j}.*RTD_cut{i,j},dt);%/...
                    %homer(RTD_cut{i,j}',dt);
                avConv_s_cut{s}(i,j) = homer(Conv_s_cut{i,j,s}.*RTD_cut{i,j}',dt);%/...
                    %homer(RTD_cut{i,j}',dt);
                avConv_mm_cut{s}(i,j) = homer(Conv_mm_cut{i,j,s}.*RTD_cut{i,j},dt);%/...
                    %homer(RTD_cut{i,j}',dt);
                
                E_cut(i,j) = homer(RTD_cut{i,j}',dt);

                %export
                VarStrArray = {'avConvL%i_pm';'avConvL%i_s';...
                    'avConvL%i_mm';'avConvL%i_exp'};
                
                VarValArray = [avConv_pm{s}(i,j),avConv_s{s}(i,j),...
                    avConv_mm{s}(i,j),avConv_exp{s}(i,j)];
                
                for p = 1:length(VarValArray)
                    varStr = sprintf(VarStrArray{p},s);
                    varVal = VarValArray(p);
                    assignin('base',varStr,varVal)
                end
                
                xavConv_pm = avConv_pm{s}(i,j);
                xavConv_s = avConv_s{s}(i,j);
                xavConv_mm = avConv_mm{s}(i,j);
                xavConv_exp = avConv_exp{s}(i,j);
                                
                tempArray(i,j) = xtemp_av;
                flowAArray(i,j) = mean(xflowA);
                k_litArray{s}(i,j) = xk_lit;
                k_expArray{s}(i,j) = mean(xk_exp(start:stop));
                taoArray(i,j) = xtao;
                tmArray(i,1) = tm(i);

            end
        end        
    end
end

% Q = row 1, T = row 3

EXPORT_XvQ         = [];
EXPORT_XvT         = [];
EXPORT_XvQcut      = [];
EXPORT_XvTcut      = [];
EXPORT_XvQideal    = [];
EXPORT_XvTideal    = [];
EXPORT_k           = [tempArray(:,3),k_expArray{1}(:,3)];

    
for s = 1:5

EXPORT_XvQ = [EXPORT_XvQ,flowAArray(:,1), avConv_exp{1}(:,1), avConv_pm{s}(:,1), avConv_s{s}(:,1),...
    avConv_mm{s}(:,1)];
EXPORT_XvT = [EXPORT_XvT,tempArray(:,3), avConv_exp{1}(:,3), avConv_pm{s}(:,3), avConv_s{s}(:,3),...
    avConv_mm{s}(:,3)];
EXPORT_XvQcut = [EXPORT_XvQcut, flowAArray(:,1), avConv_exp_cut{1}(:,1),avConv_pm{s}(:,1), avConv_s_cut{s}(:,1),...
    avConv_mm_cut{s}(:,1)];
EXPORT_XvTcut = [EXPORT_XvTcut, tempArray(:,3), avConv_exp_cut{1}(:,3), avConv_pm{s}(:,3), avConv_s_cut{s}(:,3),...
    avConv_mm_cut{s}(:,3)];
EXPORT_XvQideal = [EXPORT_XvQideal,flowAArray(:,1), avConv_exp_pm{1}(:,1), avConv_pm{s}(:,1), avConv_s_pm{s}(:,1),...
    avConv_mm_pm{s}(:,1)];
EXPORT_XvTideal = [EXPORT_XvTideal, tempArray(:,3), avConv_exp_pm{1}(:,3),avConv_pm{s}(:,3), avConv_s_pm{s}(:,3),...
    avConv_mm_pm{s}(:,3)];
EXPORT_k = [EXPORT_k k_litArray{s}(:,3)];

end

EXPORT_taoVtm = [flowAArray(:,1),taoArray(:,2),tmArray,tm_cut',tm_shift',tm_shiftcut'];
pERR = (taoArray(:,2)-tm_cut')./taoArray(:,2);


