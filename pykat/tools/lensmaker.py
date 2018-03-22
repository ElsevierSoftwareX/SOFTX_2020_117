def lensmaker(R1, R2, d, n=1.44963):
    """
    Using the Lensmaker's equation to compute the focal length of a lens with with
    thickness d, refractive index n and surface radii of curvature R1 and R2.

    Sign of R1 and R2:
    ------------------
    This function assumes that the beam travels from R1 to R2. Positive R means that
    the surface's center of curvature is further along in the direction of propagation.
    Negative R indicates that a beam that reaches the surface already has passed the
    center of curvature. Thus,
    
    R1 > 0 ---> Convex 
    R2 < 0 ---> Convex
    R1 < 0 ---> Concave
    R2 > 0 ---> Concave

    Inputs:
    -------

    R1  - Radius of curvature of the first surface in the direction of propagation [m].
    R2  - Radius of curvature of the first surface in the direction of propagation [m].
    d   - Thickness. Distance between the surfaces, measured at the optical axis [m].
    n   - Refractive index of the lens.

    Returns:
    --------
    f   - Focal length of the lens [m]. Positive for a converging lens, and negative
          for a diverging lens.  
    """
    
    return 1.0/((n-1.0)*(1.0/R1 - 1.0/R2 + (n-1.0)*d/(n*R1*R2)))

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('R1', type = float)
    parser.add_argument('R2', type = float)
    parser.add_argument('d', type = float)
    parser.add_argument('-n', '--n', type = float, required=False)

    args = vars(parser.parse_args())
    R1 = args['R1']
    R2 = args['R2']
    d = args['d']
    n = args['n']

    if n is None:
        f = lensmaker(R1, R2, d)
    else:
        f = lensmaker(R1, R2, d, n)
        

        
    
        
