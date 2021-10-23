classdef MatrixElements
    %Calculate OL,V,T between states pair.
    
    properties
        L, nPairs; Mass;
        PM, PV,parity, nPerm;
        nConf, bbox;
        data; TBP, N;
        U, Ured; %The Jacobi Coordinates tranformation matrix with N and reduced matrix in size N-1xN-1 
    end
    
    methods
        function obj = MatrixElements(data)
            %Initilize the class from the data.
            obj.L = data.BoxSize; obj.data = data;
            obj.nPairs = data.npar*(data.npar-1)/2;
            obj.nPerm = factorial(data.npar);
            obj.nConf = (2 * data.dmax + 1)^(data.npar-1);
            
            obj.Mass = diag(data.xm);
            %create TBP - Two Body Pair set of nPair vectors.
            obj.TBP = zeros(obj.nPairs, data.npar-1);
            ipair = 1; obj.N = data.npar;
            for i = 1:data.npar+1
                for j = i+1:data.npar
                    obj.TBP(ipair, i) = 1; 
                    obj.TBP(ipair, j) = -1;
                    ipair = ipair + 1;
                end
            end
            %create bbox - all combination of N particle shifts for 1 axis.
            b = zeros(1, data.npar-1) - data.dmax; 
            obj.bbox = zeros(obj.nConf, data.npar-1);
            for ibbox = 1:obj.nConf - 1
                obj.bbox(ibbox, :) = b;
                b(1) = b(1) + 1;
                for i = 1:data.npar-1
                    if b(i) == data.dmax + 1
                        b(i+1) = b(i+1) + 1;
                        b(i) = -data.dmax;
                    end
                end
            end
            obj.bbox(obj.nConf, :) = b;
            
            %create perm - all combination of N particle permutaion.
            obj.PV = perms(1:data.npar); obj.PM = zeros(data.npar, data.npar, obj.nPerm);
            for i = 1:obj.nPerm
                for j = 1:data.npar
                    obj.PM(j,obj.PV(i, j), i) = 1;
                end
                obj.parity(i) = det(obj.PM(:,:, i));
            end
                        obj.U = zeros(obj.N);
            for i = 1:obj.N-1
                for j = 1:i
                	 obj.U(i,j) = 1/sqrt(i*(i+1));
                end
                obj.U(i,i+1) = -sqrt(i/(i+1));
            end
            for i = 1:obj.N
                obj.U(obj.N,i) = 1/sqrt(obj.N);
            end
            obj.Ured = obj.U(1:obj.N-1,1:obj.N-1);
            
            for ibbox = 1:obj.nConf
                obj.bbox(ibbox,:) = (obj.Ured*obj.bbox(ibbox,:)')';
            end
        end
        
        function OL = OL(obj, state1, state2)
            %Calcualte the overlap between 2 states with eq (11) from the
            %ReseachProposalStageB by Moti
            OL = 0;
            for iperm = 1:obj.nPerm
                OLAxis = 0;
                SISPart = obj.Operator(state1, state2, 1, 2, iperm, 0);
                if SISPart == 0
                    continue;
                end
                A1_0 = state1.A(:,:);
                A2_0 = obj.PM(:,:,iperm)'*state2.A(:,:)* obj.PM(:,:,iperm);
                
                D1 = obj.U*A1_0*obj.U'; D1 = D1(1:obj.N-1,1:obj.N-1);
                D2 = obj.U*A2_0*obj.U'; D2 = D2(1:obj.N-1,1:obj.N-1);
                C = D1 + D2; CInv = C^-1;
                for iconf = 1:obj.nConf
                     v = obj.L*D1*obj.bbox(iconf,:)';
                     x4 = 0.5*v'*CInv*v;
                     x6 = -0.5*obj.L^2*obj.bbox(iconf,:)*D1*obj.bbox(iconf,:)';
                     OLAxis = OLAxis + exp(x4+x6);
                end
                OLAxis = OLAxis/sqrt(det(C));
                OL = OL + obj.parity(iperm)*SISPart*OLAxis^3;
            end
        end
        
        function E = energy(obj, state1, state2)
            %Calcualte the <psi1|H|psi2> ME between 2 states by summing 
            %eq. (12)+(13)+(14) from ReseachProposalStageB.pdf by Moti
            KinEnergy = 0; PotEnergy2B = 0; PotEnergy3B = 0;
            B = zeros(2); I = eye(2);
            for iperm = 1:obj.nPerm
                PotEnergy3BA = zeros(obj.N-2,obj.N-2,obj.N-2,3);
                PotEnergyA = zeros(obj.nPairs,obj.data.npt);
                SISPot2BPart = zeros(obj.nPairs,obj.data.npt);
                SISKinPart = obj.Operator(state1, state2, 1, 2, iperm, 0);
                OLAxis = 0; KinAxis = 0;

                A1_0 = state1.A(:,:);
                A2_0 = obj.PM(:,:,iperm)'*state2.A(:,:)* obj.PM(:,:,iperm);
                
                D1 = obj.U*A1_0*obj.U'; D1 = D1(1:obj.N-1,1:obj.N-1);
                D2 = obj.U*A2_0*obj.U'; D2 = D2(1:obj.N-1,1:obj.N-1);
                C = D1 + D2; CInv = C^-1;
                sqInvDetC = 1/sqrt(det(C));

                %%%         Calculate kinetic energy
                if SISKinPart ~= 0
                    y1 = trace(D1*CInv*D2);
                    for iconf = 1:obj.nConf
                         v = obj.L*D1*obj.bbox(iconf,:)';
                         yy = D2*CInv*v; y2 = yy'*yy;
                         x4 = 0.5*v'*CInv*v;
                         x6 = -0.5*obj.L^2*obj.bbox(iconf,:)*D1*obj.bbox(iconf,:)';
                         y3 = sqInvDetC*exp(x4+x6);
                         OLAxis = OLAxis + y3;
                         KinAxis = KinAxis + (y1-y2)*y3;
                    end
                end

                %%%         Calculate 2 body potential energy
                for ipair = 1:obj.nPairs
                    i = find(obj.TBP(ipair,:)==1); j = find(obj.TBP(ipair,:)==-1);

                    for ipt = 1:obj.data.npt
                        for iop = 1:obj.data.nop
                            SISPot2BPart(ipair,ipt) = SISPot2BPart(ipair,ipt) +...
                                obj.data.vpot(ipt,iop)*obj.Operator(state1, state2, i, j, iperm, iop-1);
                        end
                        if SISPot2BPart(ipair,ipt) ~= 0
                            GU = obj.U*obj.TBP(ipair,:)'; GU = GU(1:obj.N-1);
                            s = GU'*CInv*GU;
                            x9 = sqInvDetC*sqrt(1.0 / (2.0*obj.data.apot(ipt)*s + 1));
                            C = obj.data.apot(ipt)/(2.0*obj.data.apot(ipt)*s + 1);
                            for iconf = 1:obj.nConf
                                v = obj.L*D1*obj.bbox(iconf,:)'; CInvd = CInv*v;
                                x4 = 0.5*v'*CInv*v;
                                x5 = -0.5*obj.L^2*obj.bbox(iconf,:)*D1*obj.bbox(iconf,:)';
                                rC = GU'*CInvd;
                                Cpot = x9*exp(x4+x5);
                                for q=-obj.data.dmax:obj.data.dmax
                                    r = -obj.L*q + rC;
                                    PotEnergyA(ipair,ipt) = PotEnergyA(ipair,ipt) + Cpot*exp(-r^2*C);
                                end
                            end
                        end
                    end
                end
                %%%         Calculate 3 body potential energy
                SIS3BPotPart = obj.data.vpot3b*obj.Operator(state1, state2, 1, 2, iperm, 0);
                if SIS3BPotPart ~= 0
                    for i = 1:obj.N-2
                        for j = i+1:obj.N-1
                            for k = j+1:obj.N
                                for cyc =1:3
                                    if cyc == 1
                                        i1 = i; j1 = j; k1=k; end
                                    if cyc == 2
                                        i1 = j; j1 = k; k1=i; end
                                    if cyc == 3
                                        i1 = k; j1 = i; k1=j; end
                                    Cik = zeros(obj.N, 1); Cik(i1)=1; Cik(k1)=-1;
                                    Cjk = zeros(obj.N, 1); Cjk(j1)=1; Cjk(k1)=-1;
                                    Cik = obj.U*Cik; Cik = Cik(1:obj.N-1);
                                    Cjk = obj.U*Cjk; Cjk = Cjk(1:obj.N-1);
                                    
                                    B(1,1) = Cik'*CInv*Cik; B(1,2) = Cik'*CInv*Cjk;
                                    B(2,1) = Cjk'*CInv*Cik; B(2,2) = Cjk'*CInv*Cjk;

                                    BB = I+2*obj.data.apot3b*B;
                                    BBB = B^-1*(I-BB^-1);
                                    xx1 = sqrt(1.0 / det(BB));
                                    Temp = zeros(1,obj.nConf);
                                    for iconf = 1:obj.nConf
                                        e = zeros(2,1);
                                        v = obj.L*D1*obj.bbox(iconf,:)'; vCinv = v'*CInv;
                                        xx2 = 0.5*vCinv*v;
                                        xx3 = -0.5*obj.L^2*obj.bbox(iconf,:)*D1*obj.bbox(iconf,:)';
                                        C3f = sqInvDetC*exp(xx2+xx3)*xx1;
                                        C3ik = vCinv*Cik; C3jk = vCinv*Cjk;
                                        for q1=-obj.data.dmax:obj.data.dmax
                                            e(1)=C3ik-obj.L*q1; CC = e(1)^2*BBB(1,1);
                                            for q2=-obj.data.dmax:obj.data.dmax
                                                e(2)=C3jk-obj.L*q2;
                                                C = CC+ e(2)^2*BBB(2,2)+2*e(1)*e(2)*BBB(1,2);
                                                Temp(iconf) = Temp(iconf) + C3f*exp(-0.5*C);
                                           end
                                        end
                                    end
                                    for iconf = 1:obj.nConf
                                        PotEnergy3BA(i,j-1,k-2,cyc) = PotEnergy3BA(i,j-1,k-2,cyc)+Temp(iconf);
                                    end
                                end
                            end
                        end
                    end
                end
                %%%         Sum All Axis Kinetic Energy
                KinEnergy = KinEnergy + 3*obj.parity(iperm)*SISKinPart*KinAxis*OLAxis^2;
                %%%         Mult All Axis Potential Energy
                for ipair = 1:obj.nPairs
                    for ipt = 1:obj.data.npt
                        PotEnergy2B = PotEnergy2B + obj.parity(iperm)*SISPot2BPart(ipair,ipt)*PotEnergyA(ipair,ipt)^3;
                    end
                end
                %%%         Counter of this permutaion
                PotEnergy3BT = 0;
                for i = 1:obj.N-2
                    for j = i+1:obj.N-1
                        for k = j+1:obj.N
                            for cyc = 1:3
                                PotEnergy3BT = PotEnergy3BT + PotEnergy3BA(i,j-1,k-2,cyc)^3;
                            end
                        end
                    end
                end
                PotEnergy3B = PotEnergy3B + PotEnergy3BT*obj.parity(iperm)*SIS3BPotPart;
            end
            KinEnergy = KinEnergy*obj.data.h2m/2;
            E = KinEnergy+PotEnergy2B+PotEnergy3B;
        end
        
        function output = Operator(obj, state1, state2, ipar, jpar, perm, kind)
            x = 0; y = 0;
            tempState2 = state2;
            %%%make exchanges by the kind of the operator between the i,j particles.
            if kind == 1 || kind == 2
                for n = 1:size(state1.spin, 1)
                    tempState2.spin(n, ipar+1) = state2.spin(n, jpar+1);
                    tempState2.spin(n, jpar+1) = state2.spin(n, ipar+1);
                end
            end
            if kind == 1 || kind == 3
                for n = 1:size(state1.isospin, 1)
                    tempState2.isospin(n, ipar+1) = state2.isospin(n, jpar+1);
                    tempState2.isospin(n, jpar+1) = state2.isospin(n, ipar+1);
                    tempState2.isospin(n, 1) = -state2.isospin(n, 1);
                end
            end
            %%%sum the coeff of spin/iso states which are the same in x/y.
            for i1 = 1:size(state1.spin, 1)
                for i2 = 1:size(tempState2.spin, 1)
                    if state1.spin(i1,2:obj.N+1) == tempState2.spin(i2,obj.PV(perm, 1:obj.N)+1)
                        x = x + state1.spin(i1, 1)*tempState2.spin(i2, 1);
                    end
                end
            end
            for i1 = 1:size(state1.isospin, 1)
                for i2 = 1:size(tempState2.isospin, 1)
                    if state1.isospin(i1,2:obj.N+1) == tempState2.isospin(i2,obj.PV(perm, 1:obj.N)+1)
                        y = y + state1.isospin(i1, 1)*tempState2.isospin(i2,1);
                    end
                end
            end
            output = x*y;
        end
    end
end

