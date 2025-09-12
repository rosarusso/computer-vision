classdef TriangularMesh < handle
    %TRIANGULARMESH Triangular mesh.
    %   Use doc TriangularMesh for help.
    
    properties (SetAccess = private)
        %Mesh vertices
        verts
        %Mesh faces
        faces
        %Mesh number of vertices
        n_verts
        %Mesh number of faces
        n_faces
    end
    
    methods
        function TM = TriangularMesh(v, f)
            %TRIANGULARMESH Creates a TriangularMesh object.
            %   TM = TriangularMesh(v, f) creates the TriangularMesh object
            %   given the mesh vertices v and faces f.
            
            TM.updatemesh(v, f);
        end
        
        function updatemesh(TM, v, f)
            %UPDATEMESH Updates the TriangularMesh object.
            %
            %   TM.updatemesh(v, f) updates the TriangularMesh object TM
            %   given the new mesh vertices v and the faces f.
            
            % Check inputs
            assert(ismatrix(v) && size(v, 2) == 3, ...
                'TRIANGULARMESH:InvalidInput', ...
                'v must be a N x 3 matrix');
            
            assert(ismatrix(f) && size(f, 2) == 3, ...
                'TRIANGULARMESH:InvalidInput', ...
                'f must be a M x 3 matrix');
            
            % Get the mesh main properties
            TM.verts = v;
            TM.faces = f;
            TM.n_verts = size(v, 1);
            TM.n_faces = size(f, 1);
            
            % Check the faces indices
            assert(min(min(f)) >= 1 && max(max(f)) <= TM.n_verts, ...
                'TRIANGULARMESH:InvalidInput', ...
                'faces indices must be 1 <= idx <= n_verts');
        end
        
        function v = getverts(TM)
            %GETVERTS Get mesh vertices.
            %
            %   v = TM.getverts returns the mesh vertices.
            
            v = TM.verts;
        end
        
        function f = getfaces(TM)
            %GETFACES Get mesh faces.
            %
            %   f = TM.getfaces returns the mesh faces.
            
            f = TM.faces;
        end
        
        function n = getnverts(TM)
            %GETNVERTS Get mesh number of vertices.
            %
            %   n = getnverts returns the mesh number of vertices.
            
            n = TM.n_verts;
        end
        
        function m = getnfaces(TM)
            %GETNFACES Get mesh number of faces.
            %
            %   m = getnfaces returns the mesh number of faces.
            
            m = TM.n_faces;
        end
    end
end
