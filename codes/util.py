"""
Utility functions
"""

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    (Reference: http://stackoverflow.com/questions/9081553/python-scatter-plot-size-and-style-of-the-marker/24567352#24567352)
    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection


def peakdet(v, delta):
    """
    Peak detection by 'difference' criterion 'delta'
     Reference: https://gist.github.com/endolith/250860
    """
    import numpy as np
    
    maxtab,mintab = [],[]
    x = np.arange(len(v))
    v = np.copy(v)
    
    mn, mx = np.inf, -np.inf
    mnpos, mxpos = np.nan, np.nan
    
    lookformax = True
    
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab)

def second_largest(numbers):
    """
    Find the second largest element and its index
    """
    count = 0
    m1 = m2 = float('-inf')
    for i in range(len(numbers)):
        x = numbers[i]
        count += 1
        if x > m2:
            if x >= m1:
                idx = i
                m1, m2 = x, m1  
            else:
                m2 = x
                idx = i
    return m2,idx if count >= 2 else None

def conv_gaus(spec_list,fwhm,norm=True,
              show_norm=False,**kwargs):
    """
    Convolving 1d Gaussian to the spectra
    
    Inputs:
    - spec_list: list of numpy.1darray
      List of spectra. Presumably each entry is 
      the spectrum along each column
      
    Parameters:
    - fwhm: float
      The FWHM of the Gaussian kernel
      Note that it is related to sigma via
       FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
       
    Options
    - norm: boolean
      If True, the convolved result will be 
      normalized by Unif (conv) Gaussian
      to handle boundary effect. Default True
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    arr = np.copy(spec_list[0])    
    xx = np.arange(-arr.shape[0]/2.,arr.shape[0]/2.)
    
    ##
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    gauss = np.exp(-(xx/sigma)**2/2)
    
    ## normalization factor
    norm = np.convolve(np.ones(arr.shape[0]),gauss,mode="same",**kwargs)
    if show_norm==True:
        plt.plot(norm); plt.title('Normalization correction'); plt.show()
    
    ## convolution 
    spec_conv = []
    for i in range(len(spec_list)):
        spec_conv.append(np.convolve(spec_list[i],gauss,mode="same") / norm)
    
    return spec_conv