from PyQt5.QtWidgets import QApplication
from SliderApp.slider_fit_app import SliderFitApp
from SliderApp.slider_refl import cPlotAndFitRefl

import Reflectivity.reflmodels.reflectivity as reflectivity

import numpy as np
import sys, lmfit

class ReflGuiApp(cPlotAndFitRefl):
    def init_data(self):
        self.qmin = 1e-2
        self.qmax = np.inf

        self.x = np.linspace(0.001, 0.1, 100)

        self.p = lmfit.Parameters()
        self.p.add('dai',              1e-4,          min=0.,   max=1e-3, vary=0)
        self.p.add('I0',               1,       min=0,  max=1.5,  vary=0)
        self.p.add('Ibg',              0,            min=0.,   max=1e-4,  vary=0)
        self.p.add('roughness',        5.82017266161, min=0.0,  max=20,   vary=0)
        self.p.add('packing_density_1',0.584971550956,min=0.0,  max=1.0,  vary=0)
        self.p.add('packing_density_2',0.2,           min=0.0,  max=1.0,  vary=0)
        self.p.add('immersion_depth',    15.0834474946, min=0.,   max=100,   vary=0)
        

        # Bruker properties
        self.dlam_ov_lam = 0.01 #? 
        self.wavelength = 1.3414 # Cu-K-alpha 
        self.beamwidth = 0.2 # 200 microns
        self.samplelen = 10 # 10mm
        
        self.xsld = np.linspace(-50, 800, 450)
        self.layer_thickness = self.xsld[1] - self.xsld[0]
        self.core_radius = 93
        self.shell_thickness = 15 
        self.sig_radius = 0#0.05

        self.particle_sld = 42.146e-6
        self.oleic_sld = 8.467e-6
        self.si_sld = 20.065e-6
        
        self.get_model(self.p)

        
    def get_model(self, p):
        x0 = [0,0]
        density = [p["packing_density_1"].value, p["packing_density_2"].value]
        immersion_depth = p["immersion_depth"].value
        
        sld = reflectivity.nanoparticle_models.sld_multi_immersed_coreshell_sphere(\
            self.xsld, x0, density, self.core_radius, self.shell_thickness,
            immersion_depth, self.sig_radius,
            self.particle_sld, self.oleic_sld, self.si_sld)
        
        roughness = p["roughness"]*np.ones(len(sld))
        thickness = self.layer_thickness*np.ones(len(sld))
        
        q = self.x
        self.ymodel = reflectivity.models.parrat(q, sld, roughness,\
                            thickness)
        
        sigQ = np.sqrt((self.dlam_ov_lam * q)**2 +\
                    (4.*np.pi/self.wavelength * p["dai"])**2)
        # self.ymodel = reflectivity.math.resolution_smear(q, self.ymodel, sigQ)
        self.ymodel = p["I0"]*self.ymodel + p["Ibg"]
        
        self.ysld = reflectivity.models.roughsld_thick_layers(self.xsld,\
                                sld, roughness, thickness)
                                
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(ReflGuiApp)
    aw.show()
    app.exec_()
