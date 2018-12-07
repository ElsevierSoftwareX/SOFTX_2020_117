import numpy as np

def all_bp_detectors(_base, param='q', direction='x'):
    """Generates a string of kat code consisting of
    beam parameter detectors for every node in `_base`.

    Parameters
    ----------
    _base : :class:`.kat`
        Base object containing the optical configuration.

    param : str, optional
        Parameter of the `bp` detector to detect (defaults to 'q').

    direction : str, optional
        Plane of detection ('x' for tangential, 'y' for sagittal,
        "both" for generating `bp` detectors for both planes).

    Returns
    -------
    out : str
        A string of the kat code for the beam parameter detectors.
    """
    code = ""
    if direction == "both": direction = ['x', 'y']
    else: direction = [direction]
    for node in _base.nodes.getNodes().keys():
        for d in direction:
            code += "bp bp{n} {d} {p} {n}\n".format(
                n=node, d=d, p=param
            )
    return code

def generate_amp_detectors(node, maxtem, f=0, ad_prefix="ad"):
    """Generates a string of kat code consisting of amplitude detectors
    for each TEM up to a given maxtem.

    Parameters
    ----------
    node : str
        Name of the node at which to place the amplitude detectors.

    maxtem : int
        Maximum mode order up to which to generate amplitude detectors.

    f : double
        Frequency offset to carrier [Hz].

    ad_prefix : str
        Prefix of amplitude detector names. These will be named
        as ``[<ad_prefix>nm, ...]``.

    Returns
    -------
    out : str
        A string of the kat code for the amplitude detectors.

    author: George Smetana
    date: 04/12/2018
    """
    temsize = maxtem + 1
    code = ""
    for n in range(temsize):
        for m in range(temsize - n):
            code += "ad {prefix}{n}{m} {n} {m} {f} {node} \n".format(
                prefix=ad_prefix, n=n, m=m, f=f, node=node
            )
    return code

def get_tem_array(out, maxtem, ad_prefix='ad', power=False, order_only=True, table=False, plot=False):
    """
    Return array of TEM mode amplitudes/powers obtained from amplitude detectors at a node.

    Parameters
    ----------
    out : :class:`.KatRun`
        Output from a Finesse simulation.

    maxtem : int
        Maximum TEM order to which detectors should be generated.

    ad_prefix : str, optional
        Prefix of amplitude detector names.

    power : bool, optional
        If `True` returns the TEM powers, else returns the TEM amplitudes.

    order_only : bool, optional
        If `True`, returns quadrature sum of modes of the same TEM order,
        else returns the full 2D array.

    table : bool, optional
        If `True`, print the output using the ``tabulate`` module.

    plot : bool, optional
        If `True`, generates a plot of the output.

    Returns
    -------
    If `order_only == True` and noaxis set:
        1D array where array[k] is amp/power in kth TEM order
    If `order_only == True` and xaxis set:
        2D array where array[k, i] is amp/power in kth TEM order of ith xaxis measurement.
    If `order_only == False` and noaxis set:
        2D array where array[n, m] is amp/power in TEM_nm mode
    If `order_only == False` and xaxis set:
        3D array where array[n, m, i] is amp/power in TEM_nm mode of ith xaxis measurement.

    author: George Smetana
    date: 04/12/2018
    """
    from tabulate import tabulate
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

    temsize = maxtem + 1
    # Depth is how many measurements each detector makes, this must be set manually to 1 for special
    # case of a single measurement with noxaxis set.
    depth = 1 if out.x.shape is () else out.x.size
    # Initialise the amplitude array with all nans
    tem_amps = np.empty((temsize, temsize, depth))
    tem_amps[:] = np.nan
    for n in range(temsize):
        for m in range(temsize - n):
            tem_amps[n, m, :] = out["{prefix}{n}{m}".format(
                prefix=ad_prefix, n=n, m=m
            )]

    if depth > 1 and (table or plot):
        # If multiple measurements are made, printing a table/plotting doesn't make sense.
        raise Exception("Table printing and plotting is not supported "
                        "for multiple outputs per detector.")

    if order_only:
        orders = list(range(temsize))
        tem_powers = tem_amps*np.conj(tem_amps)

        x = np.array(orders)
        y = x[::-1]
        # Rotate array so that diagonal elements (elements of the same TEM order) are stacked in columns.
        rotated_powers = tem_powers[x, np.add.outer(x, y) - (temsize - 1), :]
        order_powers = np.nansum(rotated_powers, axis=1)
        order_vals = np.sqrt(order_powers) if not power else order_powers

        if depth == 1:
            # If only one measurement is made per detector, keeping a 2D array is redundant.
            order_vals = order_vals[:, 0]

        if table and depth == 1:
            value_name = ['Power'] if power else ['Amplitude']
            headers = ['Order'] + orders
            row = value_name + list(order_vals)
            print(tabulate([row], headers=headers, tablefmt='grid'))

        if plot and depth == 1:
            plt.figure()
            ax = plt.axes()
            ax.scatter(orders, order_vals)
            ax.set_xticks(orders)
            ax.grid(True, 'both', 'both')
            ax.set_yscale('log')
            ax.set_ylim((min(order_vals[order_vals > 0])*0.8, max(order_vals)*1.2))
            ax.set_xlabel('TEM modes order')
            ax.set_ylabel('Power' if power else 'Amplitude')
            ax.set_title('{} of Each TEM Order'.format('Power' if power else 'Amplitude'))
            plt.show()

        return order_vals
    else:
        if depth == 1:
            # If only one measurement is made per detector, keeping a 3D array is redundant.
            tem_amps = tem_amps[:, :, 0]
        tem_vals = tem_amps*np.conj(tem_amps) if power else tem_amps
        orders = list(range(temsize))

        if table and depth == 1:
            value_name = 'Power' if power else 'Amplitude'
            print('{} of All Modes - Location (n, m) for TEM_nm'.format(value_name))
            print(tabulate(tem_vals, tablefmt='grid'))

        if plot and depth == 1:
            m_list = orders * temsize
            n_list = sorted(m_list)
            plt.figure()
            ax = plt.axes()
            grid = ax.scatter(n_list, m_list, c=tem_vals.flatten(), norm=mpl.colors.LogNorm())
            plt.colorbar(grid)
            ax.set_xticks(orders)
            ax.set_yticks(orders)
            ax.grid(True, 'major', 'both')
            ax.set_xlabel('n index in $TEM_{nm}$')
            ax.set_ylabel('m index in $TEM_{nm}$')
            ax.set_title('{} of Each TEM Mode'.format('Power' if power else 'Amplitude'))
            plt.show()

        return tem_vals