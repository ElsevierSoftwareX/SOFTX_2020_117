import numpy as np
from scipy.integrate import newton_cotes

def newton_weights(x, newtonCotesOrder):
    """
    Constructs the weights for a composite Newton-Cotes rule for integration.
    These weights should be multipled by the step-size, equally spaced steps
    are assumed. The argument x should be a numpy array of the x samples points.
    
    If the order Newton-Cotes order specified does not produce a rule that fits
    into len(x), the order is decremented and that is used to fill gaps that do no fit.
    """
    # ensure we have a 1xN array
    x = np.array(x).squeeze()
    
    if newtonCotesOrder == 0:
        return np.ones(x.shape)
    
    W = newton_cotes(newtonCotesOrder, 1)[0]

    wx = np.zeros(x.shape, dtype=np.float64)    

    N = len(W)

    i = 0

    while i < len(x) - 1:
        try:
            wx[i:(i+len(W))] += W
        
            i += len(W)-1
        except ValueError as ex:
        
            if len(wx[i:(i+len(W))]) < len(W):
                newtonCotesOrder -= 1
                W = newton_cotes(newtonCotesOrder, 1)[0]
        
    return wx
    
    