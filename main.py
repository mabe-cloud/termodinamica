from BlackOil import GasPhase_Correlations
from BlackOil import OilPhase_Correlations
from barril.units import Scalar
from barril.units import ChangeScalars

# Dados do problema
T = Scalar(122, 'degF')
P = Scalar(3626, 'psia')
Pb = Scalar(5000, 'psia')
do = 0.86
dg = 0.84

T = Scalar(T.GetValue('degR'),'degR')
print(T)

print('### Gas ###')
GasPhase_Correlations(dg,P,T).output(dg,P,T)
print('')
"""
print('### Oil ###')
OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).output(do, dg,P, Pb, T)
"""