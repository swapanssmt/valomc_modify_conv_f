#ifndef _ERRORS_H_
#define _ERRORS_H_

enum mcerror {
   NO_LIGHTSOURCE,
   SIZE_MISMATCH_G,
   SIZE_MISMATCH_MUA_EX_SOL,
   SIZE_MISMATCH_MUA_EX_F,
   SIZE_MISMATCH_MUA_EM_SOL,
   SIZE_MISMATCH_MUS_EX,
   SIZE_MISMATCH_MUS_EM,
   SIZE_MISMATCH_N,
   SIZE_MISMATCH_BCN,
   MISSING_BCN,
   SIZE_MISMATCH_BCTYPE,
   SIZE_MISMATCH_LIGHT_DIRECTION,
   SIZE_MISMATCH_LIGHT_DIRECTION_TYPE,
   INCONSISTENT_H,
   INCONSISTENT_BH,
   INCONSISTENT_MESH,
   INCONSISTENT_MESH_DUPLICATE_NEIGHBORS,
   NO_GAUSSIAN_SIGMA,
   MISSING_BOUNDARY
};

// [AL]
static const char *errorstring(mcerror error) {
    switch(error) {
        case NO_LIGHTSOURCE:
            return "No lightsource";
        case SIZE_MISMATCH_G:
            return "Scattering anisotropy array is not equal in size to the number of elements";
        case SIZE_MISMATCH_MUA_EX_SOL:
            return "Absorption array of mixture at excitation wavelenght is not equal in size to the number of elements";
        case SIZE_MISMATCH_MUA_EX_F:
            return "Absorption array of fluorofore at excitation wavelenght is not equal in size to the number of elements";
        case SIZE_MISMATCH_MUA_EM_SOL:
            return "Absorption array of mixture at emission wavelenght is not equal in size to the number of elements";
        case SIZE_MISMATCH_MUS_EX:
            return "Scattering array at excitation wavelenght is not equal in size to the number of elements";
        case SIZE_MISMATCH_MUS_EM:
            return "Absorption array at emission wavelenght is not equal in size to the number of elements";
        case SIZE_MISMATCH_N:
            return "Refractive index array is not equal in size to the number of elements";
        case SIZE_MISMATCH_BCN:
            return "External refractive index array is not equal in size to the number of boundary elements";
        case MISSING_BCN:
            return "No external refractive index array given";
        case SIZE_MISMATCH_BCTYPE:
            return "The boundary condition array is not equal in size to the number of boundary elements";
        case SIZE_MISMATCH_LIGHT_DIRECTION:
            return "The light direction array is not equal in size to the number of boundary elements";        
        case SIZE_MISMATCH_LIGHT_DIRECTION_TYPE:
            return "The light direction type array is not equal in size to the number of boundary elements";        
        case INCONSISTENT_H:
            return "H matrix contain indices that are not in the coordinate matrix";                    
        case INCONSISTENT_BH:
            return "BH matrix contain indices that are not in the coordinate matrix";            
        case INCONSISTENT_MESH:
            return "Ill defined mesh - cannot find nodes. Check Mesh.";  
        case INCONSISTENT_MESH_DUPLICATE_NEIGHBORS:
            return "Ill defined mesh - some boundary elements are trapped between two elements. A boundary element must belong to a single element. Check Mesh.";  
        case MISSING_BOUNDARY:
            return "Ill defined mesh - a complete boundary cannot be constructed. Check Mesh.";                        
        case NO_GAUSSIAN_SIGMA:
            return "A gaussian light source exists but no sigma parameter was provided";                        
        default:
            return "Undefined error";
    }
}  

#endif
