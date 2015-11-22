classdef simulation < handle
    % SIMULATION Manages properties of the simulation.
    %   Specifically, SIMULATION abstracts working with the grid,
    %   so that numerical routines can be written in a simpler way, without
    %   worrying about details of the grid structure.
    properties(SetAccess = private)
        % The MSRT grid
        grid;
        
        % Number of points in x and z directions
        Nx;
        Nz;
        
        % Gravitational acceleration
        g;
        
        % Face, cell and node indices
        top_faces
        bottom_faces
        left_faces
        right_faces
        top_cells
        
        % NB! top_nodes must be sorted according to x-coordinates
        top_nodes
        
        % Initial surface height and seabed functions
        eta0
        h
    end
    
    methods
        function obj = simulation(Nx, Nz, g, eta0, h)
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
            obj.g = g;
            obj.h = h;
            obj.eta0 = eta0;
            
            [G, top, bottom, left, right] = setup_grid(eta0, h, Nx, Nz);
            obj.grid = G;
            obj.organize_indices(top, bottom, left, right);
        end
        
        function update_surface(obj, surface_shape)
            % UPDATE_SURFACE Updates the grid to match the new surface
            % shape.
            %   sim.UPDATE_SURFACE(shape) updates the grid according to the
            %   specified shape, where shape is a vector of z-values
            %   indexed such that surface_shape(top_faces(i)) is the new
            %   value for eta at the face top
            [ obj.grid, ...
                top, ...
                bottom, ...
                left, ...
                right ] ...
                = ...
                update_geometry(...
                surface_shape, ...
                obj.h, ...
                obj.Nx, ...
                obj.Nz);
            
            obj.organize_indices(top, bottom, left, right);
        end
        
        function [shape, top] = surface_shape(obj)
            % SURFACE_SHAPE Retrieve the eta values for the surface nodes.
            
            top = obj.top_nodes;
            shape = obj.grid.nodes.coords(top, 2);
        end
        
        function [grad, top] = surface_potential_gradient(obj, v)
            % SURFACE_POTENTIAL_GRADIENT Compute gradient at surface nodes
            % given half face gradient projections v.
            
            top = obj.top_nodes;
            grad = node_boundary_gradients(obj.grid, v, top);
        end
        
        function [potential, top] = surface_potential(obj, phi)
            % SURFACE_POTENTIAL Compute potential at surface nodes.
            top = obj.top_nodes;
            potential = node_boundary_potential(obj.grid, phi, top);
        end
        
        function eta = eta(obj)
            % ETA Compute surface height at surface nodes.
            eta = obj.grid.nodes.coords(obj.top_nodes, 2);
        end
        
    end
    
    methods(Access = private)
        function obj = organize_indices(obj, top, bottom, left, right)
            G = obj.grid;
            
            % Recover surface cells and nodes
            top_cell_indices = boundary_cells(G, top);
            top_node_indices = face_nodes(G, top);
            
            % Sort all indices according to their x coordinates.
            % This isn't strictly necessary but it makes it easier
            % to reason about certain things. For one, it's very convenient
            % to let eta be a vector indexed along the x coordinates of the
            % surface, so that eta(1) is the left-most node, eta(end) is
            % the right-most node.
            
            top_cells_x = G.cells.centroids(top_cell_indices, 1);
            top_nodes_x = G.nodes.coords(top_node_indices, 1);
            top_faces_x = G.faces.centroids(top, 1);
            bottom_faces_x = G.faces.centroids(bottom, 1);
            left_faces_z = G.faces.centroids(left, 2);
            right_faces_z = G.faces.centroids(right, 2);
            
            obj.top_faces = obj.sort_by_coords(top, top_faces_x);
            obj.top_cells = obj.sort_by_coords(top_cell_indices, top_cells_x);
            obj.top_nodes = obj.sort_by_coords(top_node_indices, top_nodes_x);
            obj.bottom_faces = obj.sort_by_coords(bottom, bottom_faces_x);
            obj.left_faces = obj.sort_by_coords(left, left_faces_z);
            obj.right_faces = obj.sort_by_coords(right, right_faces_z);
        end
    end
    
    methods(Static, Access = private)
        function C = sort_by_coords(C, X)
            [~, I] = sort(X);
            C = C(I);
        end
    end
    
end