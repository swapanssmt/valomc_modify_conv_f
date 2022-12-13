function vmcmedium_out = createMedium(vmcmesh, vmcmedium)
%CREATEMEDIUM Creates a vmcmedium structure
%
% USAGE:
%
%       vmcmedium = createMedium(vmcmesh);
%
%          returns
%
%                   vmcmedium.refractive_index(:) = 1;
%                   vmcmedium.scattering_coefficient_ex(:) = 0;
%                   vmcmedium.absorption_coefficient_ex_sol(:) = 0;
%                   vmcmedium.absorption_coefficient_ex_f(:) = 0;
%                   vmcmedium.absorption_coefficient_em_sol(:) = 0;
%                   vmcmedium.scattering_coefficient_em(:) = 0;
%                   vmcmedium.scattering_anisotropy(:) = 0;
%
%       vmcmedium_out = createMedium(vmcmesh, vmcmedium);
%
%           repeats the entries in vmcmedium so that the size of each array
%           is equal to the number of elements in the mesh.
%
% DESCRIPTION:
%
%       The purpose of this function is to create or resize the arrays in
%       vmcmedium so that they match the number of elements in the mesh.
%
% INPUT:
%
%       vmcmesh       -  mesh structure, contains the geometry of the
%       system
%
% OPTIONAL INPUT:
%
%       vmcmedium     -  medium structure, used to set optical
%                        coefficients
%
% OUTPUT:
%
%       vmcmedium_out -  formatted medium structure
%
% SEE ALSO:
%
%       https://inverselight.github.io/ValoMC/structures.html
%
% This function is part of ValoMC toolbox

    if(~exist('vmcmedium'))
         vmcmedium.refractive_index = 1;
         vmcmedium.scattering_coefficient_ex = 0;
         vmcmedium.absorption_coefficient_ex_sol = 0;
         vmcmedium.absorption_coefficient_ex_f = 0;
         vmcmedium.absorption_coefficient_em_sol = 0;
         vmcmedium.scattering_coefficient_em = 0;
         vmcmedium.scattering_anisotropy = 0;
    end
    vmcmedium_out = vmcmedium;

    if(~isfield(vmcmedium, 'refractive_index'))
       vmcmedium.refractive_index = 1;
    end;
    if(~isfield(vmcmedium, 'scattering_coefficient_ex'))
         vmcmedium.scattering_coefficient_ex = 0.0; 
    end;
    if(~isfield(vmcmedium, 'absorption_coefficient_ex_sol'))
         vmcmedium.absorption_coefficient_ex_sol = 0.0;
    end;
    if(~isfield(vmcmedium, 'absorption_coefficient_ex_f'))
         vmcmedium.absorption_coefficient_ex_f = 0.0;
    end;
    if(~isfield(vmcmedium, 'absorption_coefficient_em_sol'))
         vmcmedium.absorption_coefficient_em_sol = 0.0;
    end;
    if(~isfield(vmcmedium, 'scattering_coefficient_em'))
         vmcmedium.scattering_coefficient_em = 0.0; 
    end;
    if(~isfield(vmcmedium, 'scattering_anisotropy'))
         vmcmedium.scattering_anisotropy = 0.0; 
    end;
    
    % check if the optical coefficients are given as multidimensional arrays

    if(~isvector(vmcmedium.refractive_index)) 
        vmcmedium_out.nx = size(vmcmedium.refractive_index,1);
        vmcmedium_out.ny = size(vmcmedium.refractive_index,2);
        if(size(vmcmedium.refractive_index, 3)>1)
		vmcmedium_out.nz = size(vmcmedium.refractive_index,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium.scattering_coefficient_ex))
        vmcmedium_out.nx = size(vmcmedium.scattering_coefficient_ex,1);
        vmcmedium_out.ny = size(vmcmedium.scattering_coefficient_ex,2);
        
        if(size(vmcmedium.scattering_coefficient_ex, 3)>1)
		vmcmedium_out.nz = size(vmcmedium.scattering_coefficient_ex,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium.absorption_coefficient_ex_sol))
        vmcmedium_out.nx = size(vmcmedium.absorption_coefficient_ex_sol,1);
        vmcmedium_out.ny = size(vmcmedium.absorption_coefficient_ex_sol,2);
        if(size(vmcmedium.absorption_coefficient_ex_sol, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.absorption_coefficient_ex_sol,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% modified %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~isvector(vmcmedium.absorption_coefficient_ex_f))
        vmcmedium_out.nx = size(vmcmedium.absorption_coefficient_ex_f,1);
        vmcmedium_out.ny = size(vmcmedium.absorption_coefficient_ex_f,2);
        if(size(vmcmedium.absorption_coefficient_ex_f, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.absorption_coefficient_ex_f,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end

    if(~isvector(vmcmedium.absorption_coefficient_em_sol))
        vmcmedium_out.nx = size(vmcmedium.absorption_coefficient_em_sol,1);
        vmcmedium_out.ny = size(vmcmedium.absorption_coefficient_em_sol,2);
        if(size(vmcmedium.absorption_coefficient_em_sol, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.absorption_coefficient_em_sol,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end

    if(~isvector(vmcmedium.scattering_coefficient_em))
        vmcmedium_out.nx = size(vmcmedium.scattering_coefficient_em,1);
        vmcmedium_out.ny = size(vmcmedium.scattering_coefficient_em,2);
        
        if(size(vmcmedium.scattering_coefficient_em, 3)>1)
		vmcmedium_out.nz = size(vmcmedium.scattering_coefficient_em,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~isvector(vmcmedium.scattering_anisotropy))
        vmcmedium_out.nx = size(vmcmedium.scattering_anisotropy,1);
        vmcmedium_out.ny = size(vmcmedium.scattering_anisotropy,2);
        if(size(vmcmedium.scattering_anisotropy, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.scattering_anisotropy,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end

    vmcmedium_out.refractive_index = duplicateArray(vmcmedium.refractive_index(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_coefficient_ex = duplicateArray(vmcmedium.scattering_coefficient_ex(:), size(vmcmesh.H,1));
    vmcmedium_out.absorption_coefficient_ex_sol = duplicateArray(vmcmedium.absorption_coefficient_ex_sol(:), size(vmcmesh.H,1));
    vmcmedium_out.absorption_coefficient_ex_f = duplicateArray(vmcmedium.absorption_coefficient_ex_f(:), size(vmcmesh.H,1));
    vmcmedium_out.absorption_coefficient_em_sol = duplicateArray(vmcmedium.absorption_coefficient_em_sol(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_coefficient_em = duplicateArray(vmcmedium.scattering_coefficient_em(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_anisotropy = duplicateArray(vmcmedium.scattering_anisotropy(:), size(vmcmesh.H,1));
   
end


function array_out = duplicateArray(array_in, desired_size, defaultvalue)
   % convert to column vector
   if(size(array_in,2) > size(array_in,1))
       array_in = transpose(array_in);
   end

   if(size(array_in,1) > desired_size)
      tmp = array_in(1:desired_size);
   else
      tmp = repmat(array_in,ceil(desired_size/size(array_in,1)),1);
   end

   array_out = tmp(1:desired_size,:);

end
