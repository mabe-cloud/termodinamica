from BlackOil import GasPhase_Correlations
from BlackOil import OilPhase_Correlations


# Dados do problema
T = 122
P = 3626

unidades = ['psia', 'F']
Pb = 5000
unidades2= ['psia', 'psia', 'F']
do = 0.86
dg = 0.84

print('### Gas ###')
GasPhase_Correlations(dg, P, T, unidades).output()
print('')

print('### Oil ###')
OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T, Units=unidades2).output()
