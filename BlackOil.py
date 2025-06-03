import numpy as np

class set_PT_units:
    """
    Classe para conversão de unidades de temperatura e pressão.
    Fornece métodos estáticos para conversão entre diferentes escalas termométricas
    e unidades de pressão comumente usadas em engenharia e ciências.
    
    Métodos:
    ----------
        Temperature(value, from_unity, to_unity): Converte valores de temperatura
        Pressure(value, from_unity, to_unity): Converte valores de pressão
        
    Parâmetros:
    ----------
            value (float): Valor numérico da temperatura a ser convertida
            from_unity (str): Unidade de origem ('C', 'F', 'K' ou 'R')
            to_unity (str): Unidade de destino ('C', 'F', 'K' ou 'R')
            
    Retorna:
    --------
            float: Valor convertido na unidade desejada
            
        Escalas suportadas:
        - Temperatura:
            - Celsius (°C)
            - Fahrenheit (°F)
            - Kelvin (K)
            - Rankine (°R)
        - Pressão:
            - psi/psia (libra-força por polegada quadrada)
            - Pa (Pascal)
            - bar
            - atm (atmosfera padrão)
        """

    @staticmethod
    def Temperature(value, from_unity, to_unity):
        
        if from_unity == 'C':
            if to_unity == 'F':
                return (value * 9/5) + 32
            elif to_unity == 'K':
                return value + 273.15
            elif to_unity == 'R':
                return (value + 273.15) * 9/5
        
      
        elif from_unity == 'F':
            if to_unity == 'C':
                return (value - 32) * 5/9
            elif to_unity == 'K':
                return (value - 32) * 5/9 + 273.15
            elif to_unity == 'R':
                return value + 459.67
        
        
        elif from_unity == 'K':
            if to_unity == 'C':
                return value - 273.15
            elif to_unity == 'F':
                return (value - 273.15) * 9/5 + 32
            elif to_unity == 'R':
                return value * 9/5
        
        
        elif from_unity == 'R':
            if to_unity == 'C':
                return (value - 491.67) * 5/9
            elif to_unity == 'F':
                return value - 459.67
            elif to_unity == 'K':
                return value * 5/9
        
        
        if from_unity == to_unity:
            return value
        
        raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

    @staticmethod
    def Pressure(value, from_unity, to_unity):
       
        
        if to_unity == 'psi':
            if from_unity == 'Pa':
                return value / 6894.76
            elif from_unity == 'bar':
                return value * 14.5038
            elif from_unity == 'atm':
                return value * 14.6959
            elif from_unity == 'psi':
                return value
        
        
        elif from_unity == 'psi':
            if to_unity == 'Pa':
                return value * 6894.76
            elif to_unity == 'bar':
                return value / 14.5038
            elif to_unity == 'atm':
                return value / 14.6959
        
        
        elif from_unity == 'bar' and to_unity == 'atm':
            return value * 0.986923
        
       
        if from_unity == to_unity:
            return value
        
        raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")

class GasPhase_Correlations:
    """
    Classe para cálculo das propriedades termodinâmicas de gases.
    Correlações implementadas:
    - Papay (1985)
    - Lee et al. (1966)
    Retorna a massa específica, a compressibilidade e a viscosidade do gás.
    
    Parâmetros
    ----------
    dg : float
        Densidade do gás (adimensional)
    P : float
        Pressão (unidades especificadas no parâmetro units)
    T : float
        Temperatura (unidades especificadas no parâmetro units)
    units : list
        Lista com unidades de [pressão, temperatura] respectivamente
    """
    def __init__(self, dg, P, T, units):
        
        if dg == None:
            self.dg = self.rho_g()/self.rho_ar
        else:
            self.dg = dg

        a, b = units
        
        self.T = set_PT_units.Temperature(T, b, 'R')
        

    
        self.P = set_PT_units.Pressure(P, a, 'psia')

        if self.dg < 0.75:
            self.Ppc = 677 + 15.0*self.dg - 37.5*(self.dg**2)
            self.Tpc = 168 + 325*self.dg - 12.5*(self.dg**2)
        else:
            self.Ppc = 706 - 51.7*self.dg - 11.1*(self.dg**2)
            self.Tpc = 187 + 330*self.dg - 71.5*(self.dg**2)
        

        self.standardconditions()

        self.Ppr = self.P/self.Ppc
        self.Tpr = self.T/self.Tpc
        self.Mg = self.dg * 28.96 #M ar
        self.units = ['psia', 'R']

    def standardconditions(self):
        self.Psc = 14.7 #psia
        self.Tsc = 60 #F
    


    def rho_g(self):
        R = 10.73 # psia·ft³/ (lb·mol·°R)
        return (self.P * self.Mg)/(self.Z() * R * self.T)
    
    def Bg(self):
        return self.Psc/(self.Tsc) * self.Z() * (set_PT_units.Temperature(self.T,'R','F'))/self.P

   ## Papay (1985)

    def Z(self, Ppr=None, Tpr=None, dg=None):
        if Ppr is not None:
            self.Ppr = Ppr
        return 1 - (3.53 * self.Ppr) / (10 ** (0.9813 * self.Tpr)) + (0.274 * self.Ppr ** 2) / (10 ** (0.8157 * self.Tpr))

    def dZdP(self, Tpr=None, Ppr=None, h=0.5*10**(-12)):

        frwrd = self.Z(self.Ppr + h, self.Tpr)
        bckwrd = self.Z(self.Ppr - h, self.Tpr)
        return (frwrd - bckwrd) / (h)


    def mu_lee1966(self):
        
            x_v = 3.448 + (986.4 / self.T) + (0.01009 * self.Mg)
            y_v = 2.4 - (0.2 * x_v)
            k_v = ((9.379 + 0.0160 * self.Mg) * (self.T ** 1.5)) / (209.2 + (19.26 * self.Mg) + self.T)

            exp = x_v * ((self.rho_g() / 62.4) ** y_v)

            return 1e-4 * k_v * np.exp(exp)
## Dempsey (1965)

    def mu_dempsey(self, dg=None, Tpr=None, Ppr=None):

        mu_g = (1.709e-5 - 2.062e-6 * self.dg) * (set_PT_units.Temperature(self.T,'R','F')) + 8.188e-3 - 6.15e-3 * np.log10(self.dg)
        print('mu_g: ', mu_g)
        a0 = -2.4621
        a1 = 2.9705
        a2 = -2.8626e-1
        a3 = 8.0542e-3
        a4 = 2.8086
        a5 = -3.4980
        a6 = 3.6037e-1
        a7 = -1.0443e-2
        a8 = -7.9339e-1
        a9 = 1.3964
        a10 = 1.4914e-1
        a11 = 4.4102e-3
        a12 = 8.3939e-2
        a13 = -1.8641e-1
        a14 = 2.0336e-2
        a15 = -6.0958e-4

        poly = (
                a0 +
                (a1 * self.Ppr) +
                (a2 * (self.Ppr ** 2)) +
                (a3 * (self.Ppr ** 3)) +
                (self.Tpr * (a4 + (a5 * self.Ppr) + (a6 * (self.Ppr ** 2)) + (a7 * (self.Ppr ** 3)))) +
                (self.Tpr ** 2) * (a8 + (a9 * self.Ppr )+ (a10 * (self.Ppr ** 2)) + (a11 * (self.Ppr ** 3))) +
                (self.Tpr ** 3) * (a12 + (a13 * self.Ppr) + (a14 * (self.Ppr ** 2)) + (a15 * (self.Ppr ** 3)))
                 )


        mu = mu_g * np.exp(poly) / self.Tpr
        return mu

    def rho(self, P, Mg, Z, T):
        R = 0
        self.ρ = P * self.Mg / Z * R * T

    def Cpr(self, Ppr=None, dZdP=None):
        return 1/self.Ppr - (1/self.Z() * self.dZdP())
    def Cg(self):
        return self.Cpr()/self.Ppc

    def output(self):
        rho_g = GasPhase_Correlations(self.dg, self.P, self.T, self.units).rho_g()
        Cg = GasPhase_Correlations(self.dg, self.P, self.T, self.units).Cg()
        mu = GasPhase_Correlations(self.dg, self.P, self.T, self.units).mu_lee1966()

        print(f'Massa específica do gás: {rho_g:.4f} lbm/ft³')  
        print(f'Compressibilidade do gás: {Cg:.4e} psi⁻¹')  
        print(f'Viscosidade do gás: {mu:.4f} cp')
        return rho_g, Cg, mu

class OilPhase_Correlations:
    """
    Classe para cálculo das propriedades termodinâmicas do óleo.
    Calcula as viscosidades dos óleos (morto e saturado), compressibilidade e massa específica.
    
    Correlações implementadas:
    
        - Beggs e Robinson
        - Standing
    
    Parâmetros
    ----------
    do : float
        Densidade do óleo (adimensional)
        Ou None se fornecido API
    dg : float
        Densidade do gás (adimensional)
    API : float
        Grau API do óleo (opcional, se não fornecido do)
    P : float
        Pressão (unidades especificadas em Units)
    Pb : float
        Pressão de bolha (unidades especificadas em Units)
    T : float
        Temperatura (unidades especificadas em Units)
    Units : list
        Lista com unidades de [P, Pb, T] (ex: ['psia', 'psia', 'F'])

    """
    def __init__(self, do, dg, API, P, Pb, T, Units):

        a, b, c = Units
        
        self.P = set_PT_units.Pressure(P, a, 'psia')
        self.T = set_PT_units.Temperature(T, c, 'F')

        self.dg = dg

        if do == None:
            self.API = API
            self.do = 141.5/(API + 131.5)
        else:
            self.do = do
            self.API = (141.5/do) - 131.5

        if Pb == None:
            self.Pb = self.bubblepressure()
        else:
            self.Pb = set_PT_units.Pressure(Pb, b, 'psia')
            self.Rs()

        self.units = ['psia', 'psia', 'F']
    
        
    ## Standing
    def bubblepressure(self):

        a = 0.00091 * self.T - 0.0125 * self.API

        return 18.2 * ((self.Rs())/self.dg) * (10**a) - 1.4


    def Rs(self):

        if self.P > self.Pb:
            a = 0.0125*self.API - 0.00091*self.T
            return self.dg * (((self.Pb/18.2) + 1.4) * 10**a)**(1/0.83)
        else:

            a = 0.0125*self.API - 0.00091*self.T
            return self.dg * (((self.P/18.2) + 1.4) * 10**a)**(1/0.83)

    def dRsdP(self):
        a = 0.0125*self.API - 0.00091*self.T
        return self.dg * (((1/18.2) + 1.4) * 10**a)**(1/0.83)
    
    def dBodP(self):
        return 0.00012*(self.dRsdP() * (1.2)*(self.dg/self.do)**(0.5) + 1.25*self.T)**(0.2)
    
    def Co(self):
        if self.P >= self.Pb:
            a = self.rho_ob() + 0.004347 * (self.P - self.Pb) - 79.1
            b = 0.0007141 * (self.P - self.Pb) - 12.938
            return 1e-6 * np.exp(a/b)
        else:
            u = ['psia', 'R']
            return -1/self.Bo() * self.dBodP() + GasPhase_Correlations(self.dg, self.P, set_PT_units.Temperature(self.T, 'F', 'R'), u).Bg()/self.Bo() * self.dRsdP()
        
    def Bo(self):
        if self.P > self.Pb:
            return self.Bob * np.exp(-self.Co() * (self.P - self.Pb))

        else:
            return 0.9759 + 0.00012*(self.Rs() * (self.dg/self.do)**(0.5) + 1.25*self.T)**(1.2)
            
    # Fim das correlações Standing
    def rho_o(self):
        if self.P > self.Pb:
            return self.rho_ob() * np.exp(self.Co()*(self.P-self.Pb))

        else:
            return (62.4 * self.do + 0.0136 * self.Rs() * self.dg)/self.Bo()
        
    def rho_ob(self):
        return Exception('Não implementado')
        
    # Correlação de Beggs e Robinson
    def mu_od(self):
        A = 10 ** (3.0324 - 0.02023 * self.API)
        return (10 ** (A * self.T ** (-1.163))) - 1
    
    def mu_ob(self):
        a = 10.715 * (self.Rs() + 100)**(-0.515)
        b = 5.44 * (self.Rs() + 150)**(-0.338)
        return a * self.mu_od()**b

    def output(self):
        
        rho_o = OilPhase_Correlations(do=self.do, dg=self.dg, API=None, P=self.P, Pb=self.Pb, T=self.T, Units = self.units).rho_o()
        Co = OilPhase_Correlations(do=self.do, dg=self.dg, API=None, P=self.P, Pb=self.Pb, T=self.T, Units = self.units).Co()
        mu_od = OilPhase_Correlations(do=self.do, dg=self.dg, API=None, P=self.P, Pb=self.Pb, T=self.T, Units = self.units).mu_od()
        mu_ob = OilPhase_Correlations(do=self.do, dg=self.dg, API=None, P=self.P, Pb=self.Pb, T=self.T, Units = self.units).mu_ob()

    
        print(f'Massa específica do óleo: {rho_o:.2f} lbm/ft³')  
        print(f'Compressibilidade do óleo: {Co:.4e} psi⁻¹')  
        print(f'Viscosidade do óleo morto: {mu_od:.3f} cp')  
        print(f'Viscosidade do óleo saturado: {mu_ob:.3f} cp')
        return rho_o, Co, mu_od, mu_ob
