%% N body penomelogical gauusian potential
%% section 1 : get parameters from txt
clearvars; close all;
fileName = "He4";
data = InputData(fileName + ".inp");
RandG = RandStream('mt19937ar','Seed',abs(2*data.irand));
ME = MatrixElements(data);
SVM = SVM(RandG, data, ME);
N = data.npar;
SVM.initilize(fileName + ".basis");
%% section 2 - SVM
for itr = SVM.StatesCounter:data.mnb
    tic
    SVM.addState();
    SVM.SaveStatesToFile(fileName + ".basis");
    msg = "itr = " + num2str(itr+1) + "    E = " + num2str(SVM.getEnergy(),'%.8f')
    toc
end

%% section3 - analyze A,B,s dis
AList = zeros(1,3*SVM.StatesCounter);
for i=1:SVM.StatesCounter
    AList(3*i-2) = SVM.States(i).A(1,2);
    AList(3*i-1) = SVM.States(i).A(1,3);
    AList(3*i) = SVM.States(i).A(2,3);
end
AdList = sqrt(-2./AList);

histogram(AdList,8)
title('AdList')
grid on; box on;
xlabel('a_{ME} with amax = 3[fm]', 'fontsize', 16); 
ylabel('counter', 'fontsize', 16); 



%% add one - save He4
basisFileName = "He4.basis";
fid = fopen(basisFileName,'rt') ;
S = textscan(fid,'%s','Delimiter','\n'); 
BasisCell = S{1};
A = zeros(880,4,4);
for itr = 1:880
    lineA = 2*itr; Dtxt = string(BasisCell(lineA));
    AME = textscan(Dtxt,'%s','Delimiter',' '); AME = AME{1};
    for i = 1:4
        for j = 1:4
            A(itr, i,j) = str2num(string(AME(3*(i-1)+j)));
        end
    end
end

for i=1:SVM.StatesCounter
    AList(6*i-5) = SVM.States(i).A(1,2);
    AList(6*i-4) = SVM.States(i).A(1,3);
    AList(6*i-3) = SVM.States(i).A(2,3);
    AList(6*i-2) = SVM.States(i).A(1,2);
    AList(6*i-1) = SVM.States(i).A(1,3);
    AList(6*i) = SVM.States(i).A(2,3);
end
AdList = sqrt(-2./AList);

histogram(AdList,8)
title('AdList')
grid on; box on;
xlabel('a_{ME} with amax = 3[fm]', 'fontsize', 16); 
ylabel('counter', 'fontsize', 16); 


%% section 4: Calculate Spin Isospin WF's
%% Calculation for S=1, T=0:
SUM = 0;
State = SVM.States(1);
for iP = 1:4
    for jP = iP+1:4
        for iperm = 1:ME.nPerm
            SUM = SUM + ME.Operator(State, State, iP, jP, iperm, 0);
            SUM = SUM - ME.Operator(State, State, iP, jP, iperm, 1);
            SUM = SUM + ME.Operator(State, State, iP, jP, iperm, 2);
            SUM = SUM - ME.Operator(State, State, iP, jP, iperm, 3);
        end
    end
end
SUM

%% Calculation for S=0, T=1:
SUM = 0;
State = SVM.States(1);
for iP = 1:4
    for jP = iP+1:4
        for iperm = 1:ME.nPerm
            SUM = SUM + ME.Operator(State, State, iP, jP, iperm, 0);
            SUM = SUM - ME.Operator(State, State, iP, jP, iperm, 1);
            SUM = SUM - ME.Operator(State, State, iP, jP, iperm, 2);
            SUM = SUM + ME.Operator(State, State, iP, jP, iperm, 3);
        end
    end
end
SUM