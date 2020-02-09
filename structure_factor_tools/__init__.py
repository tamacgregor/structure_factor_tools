import sys

__version__ = '0.0.1'


# python version
_p3 = False
if sys.version_info > (3, 0):
    _p3 = True
del(sys)

# old hyperspy wants API 'QDate' set to version 2
if not _p3:
    try:
        import sip
        sip.setapi('QDate', 2)
        sip.setapi('QDateTime', 2)
        sip.setapi('QString', 2)
        sip.setapi('QTextStream', 2)
        sip.setapi('QTime', 2)
        sip.setapi('QUrl', 2)
        sip.setapi('QVariant', 2)
        del(sip)
    except:
        pass

__all__ = [
    'structure_factor_tools_calculate',

    ]


# To get sub-modules
for x in __all__:
    exec('from . import %s' %(x))
del(x)

# Import classes
from .structure_factor_tools_calculate import StructureFactorSimulation
del(StructureFactorSimulation)
