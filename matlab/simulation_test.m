classdef simulation_test < matlab.unittest.TestCase
    
    properties
        sim
        surface_x
    end
    
    methods(TestMethodSetup)
        function setupSim(test)
            L = 5;
            D = 5;
            nx = 50;
            nz = 50;
            eta = @(x) zeros(size(x));
            h = @(x) zeros(size(x));
            
            s = simulation(L, D, nx, nz, eta, h);
            top = s.surface_nodes;
            test.surface_x = s.grid.nodes.coords(top, 1);
            test.sim = s;
        end
    end
    
    methods(TestMethodTeardown)
        
    end
    
    methods (Test)
        
        function update_surface_yields_expected(test)
            S = test.sim;
            X = test.surface_x;
            
            new_surface = sin(5 * pi * X);
            S.update_surface(new_surface);
            test.assertEqual(S.surface_shape(), new_surface, 'AbsTol', 1e-12);
        end
    end
    
end

