"""
Prints an array, with columns properly aligned, to the Python's console.
It is very useful during the code development.

Individual entries are possible to print with default format (hard-coded), or
by providing additional string as argument with desired format.

Number of arguments can be one (for default format), or two, if format is
specified with the function call.
"""

# =============================================================================
def print_array(*args):
# -----------------------------------------------------------------------------
    """
    Args: Number of input arguments can be one or two, depending if one
          wants to print the array with entries in default format (one
          argument), or with specified format (two arguments).

      One argument provided (for defualt format)
        args[0]: Array to be printed.

      Five arguments provided (for non-uniform meshes)
        args[0]: Same as above.
        args[1]: String with desired format.

    Returns:
      none!
    """

    if len((args)) == 1:
        a      = args[0]
        format = "%7.2f"
    else:
        a      = args[0]
        format = args[1]

    print('Array ['+('%d' %a.shape[0])+']')
    rows = a.shape[0]
    for i in range(0,rows):
        print(format %a[i])

    return  # end of function
