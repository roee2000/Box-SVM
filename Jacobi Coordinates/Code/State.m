classdef State
    %State of N particles with the gauusian structure of 3 axis and |S,T>
    
    properties
        A;
        spin, isospin;
    end
    
    methods
        function obj = State(A, spin, isospin)
            %Copy state properties
            obj.A = A;
            obj.spin = spin; obj.isospin=isospin;
        end
    end
end

