"""
Elegantly prints the calling function name on the Python console.
"""

# =============================================================================
def at(f_name):
# -----------------------------------------------------------------------------
    """
    Args:
      f_name: Function name.

    Returns:
      None!
    """

#    # Version 0
#    print("+" + "-" * (len(f_name) + 8) + "+")
#    print("| From: " + f_name + " |")
#    print("+" + "-" * (len(f_name) + 8) + "+")

    # Version 1
    print("=---> " + f_name)
    
    return  # end of function
