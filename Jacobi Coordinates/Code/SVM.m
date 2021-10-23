classdef SVM < handle
    %Manage the SVM procedure
    
    properties
        data
        RandG, ME
        N, seed, mm0, kk0, mnb;
        bmin, bmax,amin,amax;
		iBoxInf;  basisFileName;
        Norm, H;
        States, StatesCounter;
        Ek, dE;
        nSConf
        U;
    end
    
    methods
        function obj = SVM(RandG, data, ME)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Ek = 0; obj.dE = 0;
            obj.RandG = RandG; obj.ME = ME;
            obj.data = data;
            obj.N = data.npar; obj.StatesCounter = 0;
            obj.amin = data.amin; obj.amax = data.amax;
            obj.mm0 = data.mm0; obj.kk0 = data.kk0;
            obj.mnb = data.mnb;
            obj.nSConf = 2^obj.N;
        end
        
        function obj = initilize(obj,basisFileName)
            %Start the states basis, if there is something in the basis
            %file - use it. Otherwise add first state.
            obj.basisFileName = basisFileName; 
            if isfile(basisFileName)
                fid = fopen(basisFileName,'rt') ;
                S = textscan(fid,'%s','Delimiter','\n'); 
                if isempty(S{1})
                    basisSize = 0;
                else
                    basisSize=str2num(string(S{1}(1)));
                end
                if basisSize > 0
                    obj.getExistingBasis(S{1});
                end
            end
        end
        
        function obj = getExistingBasis(obj,BasisCell)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here basisFileName
            A = zeros(obj.N);
            for itr = 1:str2num(string(BasisCell(1)))
                tic
                obj.StatesCounter = itr;
                lineA = 2*itr; Dtxt = string(BasisCell(lineA));
                AME = textscan(Dtxt,'%s','Delimiter',' '); AME = AME{1};
                for i = 1:obj.N
                    for j = 1:obj.N
                        A(i,j) = str2num(string(AME(obj.N*(i-1)+j)));
                    end
                end
                if itr == 1
                    obj.States = State(A, obj.data.spin, obj.data.isospin);
                end
                obj.States(itr) = State(A, obj.data.spin, obj.data.isospin); 
                obj.updateNorm(obj.States(itr)); obj.updateH(obj.States(itr));
                msg = "itr = " + num2str(itr) + "    E = " + num2str(obj.getEnergy(),'%.8f')
                toc
            end
        end
        
        function updateNorm(obj, state)
            %calculate ME for Norm with the new state
            indx = obj.StatesCounter;
            obj.Norm(indx, indx) = obj.ME.OL(state, state);
            for i = 1:indx-1
                obj.Norm(i, indx) = obj.ME.OL(obj.States(i), state);
                obj.Norm(indx, i) = obj.Norm(i, indx);
            end
        end
        
        function obj = updateH(obj, state)
            %calculate ME for Norm with the new state
            indx = obj.StatesCounter;
            obj.H(indx, indx) = obj.ME.energy(state, state);
            for i = 1:indx-1
                obj.H(i, indx) = obj.ME.energy(obj.States(i), state);
                obj.H(indx, i) = obj.H(i, indx);
            end    
        end
        
        function e = getEnergy(obj)
            %save the lowest enery for k states in Ek(k). every 5 state
            %print long msg
            E = eigs(obj.H, obj.Norm, obj.StatesCounter);
            %Update Ek and dE for the selection process
            obj.Ek(obj.StatesCounter) = min(E);
            if obj.StatesCounter>3
                obj.dE = (obj.Ek(obj.StatesCounter-3)-obj.Ek(obj.StatesCounter))/3;
            else
                obj.dE = 1;
            end
            e = min(E);
        end
        
        function obj = SaveStatesToFile(obj,basisFileName)
            %Save last state to .basis file
            fileID = fopen(basisFileName,'w');
            fprintf(fileID, num2str(obj.StatesCounter));
            
            for c = 1:obj.StatesCounter
                state = obj.States(c);
                fprintf(fileID, "\n");
                %print A
                for i = 1:obj.N
                    for j = 1:obj.N
                        fprintf(fileID, num2str(state.A(i, j),'%20.12f'));
                        fprintf(fileID, " ");
                    end
                end
                fprintf(fileID, "\n");
            end
            fclose(fileID);
        end
        
        function [newState, A] = getRandomState(obj)
            %Add a new state and update the Norm and H matrices
            %Create A Matrix
            dij = zeros(obj.N);
            for i = 1:obj.N
                for j = i+1:obj.N
                    dij(i,j) = obj.RandG.rand()*(obj.data.amax-obj.data.amin)+obj.data.amin;
                    dij(j,i)=dij(i,j);
                end
            end
            A = obj.getA(dij);            
            newState = State(A, obj.data.spin, obj.data.isospin);
            
            obj.updateNorm(newState);
            if obj.validNorm() == 0 
                [newState, Ddij] = obj.getRandomState();
            end
        end
        
        function obj = addState(obj)
            obj.data.dmax = obj.data.dmaxLow;
            obj.ME = MatrixElements(obj.data);
            
            obj.StatesCounter = obj.StatesCounter + 1;
            [newStateTemp,Ddij] = obj.getRandomState();
            obj.updateNorm(newStateTemp); 
            obj.updateH(newStateTemp);
            DT = zeros(obj.N, obj.N, obj.data.mm0+1);
            
            E(1) = obj.getEnergy();
            for k = 1:obj.data.kk0
                for c1 = 1:obj.N
                    for c2 = c1+1:obj.N
                        DT(:,:,1) = Ddij;
                        for m = 2:obj.data.mm0+1
                            DT(:,:,m) = Ddij;
                            DT(c1,c2,m) = obj.RandG.rand()*(obj.data.amax-obj.data.amin)+obj.data.amin;
                            DT(c2,c1,m)=DT(c1,c2,m);
                            A = obj.getA(DT(:,:,m));
                            newStateTemp(m) = newStateTemp(1); 
                            newStateTemp(m).A = A;
                            obj.updateNorm(newStateTemp(m));
                            if obj.validNorm()
                                obj.updateH(newStateTemp(m));
                                E(m) = obj.getEnergy();
                            else
                                E(m) = inf;
                            end
                        end
                        newStateTemp(1) = newStateTemp(E==min(E));
                        Ddij = DT(:,:,E==min(E));
                        E(1) = min(E);
                    end
                end
            end
            if obj.StatesCounter == 1
                obj.States = newStateTemp(1);
            else
                obj.States(obj.StatesCounter) = newStateTemp(1);
            end
            obj.data.dmax = obj.ME.data.dmaxHigh;
            obj.ME = MatrixElements(obj.data);
            
            obj.updateNorm(newStateTemp(1)); obj.updateH(newStateTemp(1));
        end
        
        function A = getA(obj, dij)
            A = zeros(obj.N);
            for i = 1:obj.N
                for j = i:obj.N
                    if i == j
                        for k = 1:obj.N
                            if i ~= k
                                A(i, j) = A(i, j) + 2 * dij(i, k)^-2;
                            end
                        end
                    else
                        A(i, j) = -2 * dij(i, j)^-2;
                        A(j, i) = A(i, j);
                    end
                end
            end
        end
        
        function output = validNorm(obj)
            %check if the Norm matrix is valid - contains no zeros eigenvalues.
            output = true; itr=obj.StatesCounter;
            for i=1:obj.StatesCounter-1
                vdotv = obj.Norm(itr, i) / sqrt(obj.Norm(i, i) * obj.Norm(itr, itr));
                if vdotv > 0.99
                    output = false;
                end
            end
            if obj.Norm(itr, itr) < 1e-7
                output = false;
            end
            mini = min(eigs(obj.Norm,obj.StatesCounter));
            maxi = max(eigs(obj.Norm,obj.StatesCounter));
            if mini < 2e-12
                output = false;
            end
            if mini/maxi < 2e-13
                output = false;
            end
        end
    end
end

