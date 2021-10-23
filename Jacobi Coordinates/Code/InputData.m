classdef InputData
    %This object hold the data from the input file (for exanple He3.inp)    
    properties
        fileName; %example ("He3.inp")
        npar; %number of particles
        xm; %masses of the particles in atomic units
        nisc; cisc; iso; %number of isospin configuration and the configuration.
        %2 particle singlet structure for example - 
        %cisc(1)=1  iso(1,1)=1  iso(2,1)=2 
        %cisc(2)=-1 iso(1,2)=2  iso(2,2)=1
        nspc; cspc; isp; %same for spin
        h2m %usually heavy proton with h2m = 23.83
        irand; ico; %ibf
        mm0; kk0; mnb; amin; amax; %SVM constants
        BoxSize; dmax; dmaxLow; dmaxHigh; %BoxSize in fm, amount of shift, precision digit
        npt; nop; %Number of potential and operators
        vpot; apot; vpot3b; apot3b;%2-3 body potential constants
        
        spin, isospin;
    end
    
    methods
        function obj = InputData(fileName)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.fileName = fileName; 
            counter = 0;
            fid = fopen(fileName,'rt') ;
            S = textscan(fid,'%s','Delimiter','\n');
            file=string(S{1}); txt = "";
            for line = 1:size(file, 1)
                lineCharArray = char(file(line, 1));
                if isempty(lineCharArray) || lineCharArray(1) == '!'
                    counter = counter + 1;
                else
                    txt = txt + file(line) + " ";
                end
            end
            
            txtChar = char(txt);
            obj.npar = obj.getParameter(txtChar,'npar');
            obj.nisc = obj.getParameter(txtChar,'nisc');
            obj.nspc = obj.getParameter(txtChar,'nspc');
            obj.h2m = obj.getParameter(txtChar,'h2m');
            obj.irand = obj.getParameter(txtChar,'irand');
            obj.ico = obj.getParameter(txtChar,'ico');
            obj.mm0 = obj.getParameter(txtChar,'mm0');
            obj.kk0 = obj.getParameter(txtChar,'kk0');
            obj.mnb = obj.getParameter(txtChar,'mnb');
            obj.amin = obj.getParameter(txtChar,'amin');
            obj.amax = obj.getParameter(txtChar,'amax');
            obj.BoxSize = obj.getParameter(txtChar,'BoxSize');
            obj.dmaxLow = obj.getParameter(txtChar,'dmaxLow');
            obj.dmaxHigh = obj.getParameter(txtChar,'dmaxHigh');
            obj.dmax = obj.dmaxHigh;
            obj.npt = obj.getParameter(txtChar,'npt');
            obj.nop = obj.getParameter(txtChar,'nop');
            obj.vpot3b = obj.getParameter(txtChar,'vpot3b');
            obj.apot3b = obj.getParameter(txtChar,'apot3b');
            
            obj.xm = obj.getVectorParameter(txtChar,'xm');
            obj.cisc = obj.getVectorParameter(txtChar,'cisc');
            obj.cspc = obj.getVectorParameter(txtChar,'cspc');
            obj.apot = obj.getVectorParameter(txtChar,'apot');
            
            obj.iso = obj.getMatrixParameter(txtChar,'iso(');
            obj.isp = obj.getMatrixParameter(txtChar,'isp(');
            obj.spin = zeros(obj.nspc, obj.npar + 1);
            for i = 1:obj.nspc
                obj.spin(i, 1) = obj.cspc(i);
                for j = 1:obj.npar
                    obj.spin(i, j+1) = obj.isp(j,i);
                end
            end
            for i = 1:obj.nisc
                obj.isospin(i, 1) = obj.cisc(i);
                for j = 1:obj.npar
                    obj.isospin(i, j+1) = obj.iso(j,i);
                end
            end
            obj.vpot = obj.getMatrixParameter(txtChar,'vpot(');
        end
        
        function value = getParameter(obj, txtChar, parameter)
            seqStart = strfind(txtChar, parameter);
            valueIdxStart = seqStart + length(parameter)+1;
            valueIdxEnd = valueIdxStart;
            while valueIdxEnd ~= length(txtChar) && txtChar(valueIdxEnd+1) ~= ' '
                valueIdxEnd = valueIdxEnd+1;
            end
            value = str2num(txtChar(valueIdxStart:valueIdxEnd));
        end
        
        function vector = getVectorParameter(obj, txtChar, parameter)
            seqStart = strfind(txtChar, parameter);
            for counter = 1:length(seqStart)
                valueIdxStart = seqStart(counter) + length(parameter)+4;
                valueIdxEnd = valueIdxStart;
                while (valueIdxEnd < length(txtChar) && txtChar(valueIdxEnd+1) ~= ' ')
                    valueIdxEnd = valueIdxEnd+1;
                end
                vector(str2num(txtChar(valueIdxStart-3))) = str2num(txtChar(valueIdxStart:valueIdxEnd));
            end
        end
        
        function matrix = getMatrixParameter(obj, txtChar, parameter)
            seqStart = strfind(txtChar, parameter);
            for counter = 1:length(seqStart)
                valueIdxStart = seqStart(counter) + length(parameter)+5;
                valueIdxEnd = valueIdxStart;
                while (valueIdxEnd < length(txtChar) && txtChar(valueIdxEnd+1) ~= ' ')
                    valueIdxEnd = valueIdxEnd+1;
                end
                i = str2num(txtChar(valueIdxStart-5)); j = str2num(txtChar(valueIdxStart-3));
                matrix(i, j) = str2num(txtChar(valueIdxStart:valueIdxEnd));
            end
        end
    end
end

