#data class

import os.path
import pylab

class Data():
    
    def __init__(self,base,**kwds):
        self.__dict__.update(kwds)
        self.base = base
        #the default name is the base, removing the path 
        self.label = self.base.split('/')[-1]
        self.output = {}
        self.update()
        
    #shortcut
    def o(self,name):
        return self.output[name]
    
    def __str__(self):
        return " model.label =" + self.label + "\n model.base = "+self.base
    
    #add here more file extensions for the output
    output_names_list = ['background','pk','cl','cl_lensed',
                         'perturbations_k0_s','perturbations_k0_t',
                         'perturbations_k1_s','perturbations_k1_t',
                         'perturbations_k2_s','perturbations_k2_t']    
    
    def update(self):
        self.update_output()
    
    def update_output(self):
        for name in self.output_names_list:
            filename = self.base+name+".dat"
            if os.path.isfile(filename):
                self.output[name] = pylab.loadtxt(filename,
                                                  comments='#')
        if len(self.output)==0:
            print "No output found in "+self.base
        #if os.path.isfile(self.base+"background.dat"):
            #self.output['background'] = pylab.loadtxt(self.base+"background.dat",
                                                      #comments='#')