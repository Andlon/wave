classdef simulation < handle
    % SIMULATION Manages properties of the simulation.
    %   Specifically, SIMULATION abstracts working with the grid,
    %   so that numerical routines can be written in a simpler way, without
    %   worrying about details of the grid structure.
    
    properties
        % Gravitational acceleration
        g = 9.81
    end
    
    properties(SetAccess = private)
        % The MSRT grid
        grid
        
        % Physical dimension in x-direction
        L = 1
        
        % Depth of water
        D = 1
        
        % Number of points in x and z directions
        Nx = 50
        Nz = 50
        
        % Surface face and node indices
        surface_faces
        surface_nodes
        
        % Initial surface height and seabed functions
        eta0
        h
    end
    
    methods
        function obj = simulation(L, D, Nx, Nz, eta0, h)
            % SIMULATION Sets up a new wave simulation.
            %   sim = SIMULATION(Nx, Nz, g, eta0, h) instantiates a
            %   simulation object with a grid of size (Nx) x (Nz), with
            %   gravitational acceleration g, initial surface height eta0
            %   and seabed shape h. Note that eta0 and h are functions,
            %   and in particular BOTH functions should be defined relative
            %   to z = 0.
            
            require mimetic;
            obj.Nx = Nx;
            obj.Nz = Nz;
            obj.h = h;
            obj.eta0 = eta0;
            
            [G, faces, nodes] = setup_grid(L, D, eta0, h, Nx, Nz);
            obj.grid = G;
            obj.surface_faces = faces;
            obj.surface_nodes = nodes;
            obj.organize_indices();
        end
        
        function update_surface(obj, surface_shape)
            % UPDATE_SURFACE Updates the grid to match the new surface
            % shape.
            %   sim.UPDATE_SURFACE(shape) updates the grid according to the
            %   specified shape, where shape is a vector of z-values
            %   indexed such that surface_shape(top_faces(i)) is the new
            %   value for eta at the face top
            % [G, surface_faces, surface_nodes] = update_geometry(obj.L, obj.D, surface_shape, obj.h, obj.Nx, obj.Nz);
            %                 surface_shape, ...
            %                 obj.h, ...
            %                 obj.Nx, ...
            %                 obj.Nz);
            
            dx = obj.L / obj.Nx;
            nodefunc = @(x) surface_shape(round(1 + x / dx));
            [G, faces, nodes] = setup_grid(obj.L, obj.D, nodefunc, ...
                obj.h, obj.Nx, obj.Nz);
            
            obj.grid = G;
            obj.surface_faces = faces;
            obj.surface_nodes = nodes;
            
            obj.organize_indices();
        end
        
        function [shape, top] = surface_shape(obj)
            % SURFACE_SHAPE Retrieve the eta values for the surface nodes.
            
            top = obj.surface_nodes;
            shape = obj.grid.nodes.coords(top, 2);
        end
        
        function [grad, top] = surface_potential_gradient(obj, v)
            % SURFACE_POTENTIAL_GRADIENT Compute gradient at surface nodes
            % given half face gradient projections v.
            
            top = obj.surface_nodes;
            grad = node_boundary_gradients(obj.grid, v, top);
        end
        
        function [potential, top] = surface_potential(obj, phi)
            % SURFACE_POTENTIAL Compute potential at surface nodes.
            top = obj.surface_nodes;
            potential = node_boundary_potential(obj.grid, phi, top);
        end
        
        function eta = eta(obj)
            % ETA Compute surface height at surface nodes.
            eta = obj.grid.nodes.coords(obj.surface_nodes, 2);
        end
        
    end
    
    methods(Access = private)
        function obj = organize_indices(obj)
            G = obj.grid;
            faces = obj.surface_faces;
            nodes = obj.surface_nodes;
            nodes_x = G.nodes.coords(nodes, 1);
            faces_x = G.faces.centroids(faces, 1);
            obj.surface_faces = obj.sort_by_coords(faces, faces_x);
            obj.surface_nodes = obj.sort_by_coords(nodes, nodes_x);
        end
    end
    
    methods(Static, Access = private)
        function C = sort_by_coords(C, X)
            [~, I] = sort(X);
            C = C(I);
        end
    end
    
end