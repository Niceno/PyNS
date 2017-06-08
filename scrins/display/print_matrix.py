"""
Prints two-dimensional matrix, with columns properly aligned, to the Python's
console.  It is very useful during the code development.

Individual entries are possible to print with default format (hard-coded), or
by providing additional string as argument with desired format.

Number of arguments can be one (for default format), or two, if format is
specified with the function call.
"""

# =============================================================================


def print_matrix(*args):
    # -----------------------------------------------------------------------------
    """
    Args: Number of input arguments can be one or two, depending if one wants
          to print the matrix with entries in default format (one argument),
          or with specified format (two arguments).

      One argument provided (for defualt format)
        args[0]: Matrix to be printed.

      Five arguments provided (for non-uniform meshes)
        args[0]: Same as above.
        args[1]: String with desired format.

    Returns:
      none!
    """

    if len((args)) == 1:
        a = args[0]
        format = "%7.2f"
    else:
        a = args[0]
        format = args[1]

    print('Matrix[' + ('%d' % a.shape[0]) + '][' + ('%d' % a.shape[1]) + ']')
    rows = a.shape[0]
    cols = a.shape[1]

    for i in range(0, rows):
        for j in range(0, cols):
            print((format % a[i, j]), end='')
        print('')
    print('')

    return  # end of function
