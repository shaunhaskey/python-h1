

"""
Some usefull and used error classes for H1 python packages/routines


"""


class Grid3DError(Exception):
    """
    Unhandled error of grid3D
    """
    pass



class UnknownCoordSysError(Exception):
    """
    Unhandled error of grid3D
    """
    pass


class InputDataError(Exception):
    """
    Unhandled error caused by input data
    """
    pass


class GeometryError(Exception):
    """
    Unhandled geometry error 
    """
    pass


class EquilibriumError(Exception):
    """
    Unknown equilibrium type
    """
    pass

