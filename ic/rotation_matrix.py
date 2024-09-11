
import sys
import numpy as np
import configparser

class RotationMatrix():
    def __init__(self, config_file, args):
        config     = configparser.ConfigParser()
        config.read(config_file)

        self.ngc_sgc                 = args.ngc_sgc

        if self.ngc_sgc is None:
            self.ngc_sgc       =  config.get('sim', 'ngc_sgc')

        if self.ngc_sgc == "NGC":
            self.rotation_matrix = self.ngc_matrix()
        elif self.ngc_sgc == "SGC":
            self.rotation_matrix = self.sgc_matrix()
        else:
            self.rotation_matrix = None
            print("ERROR: wrong chosen galactic cap")
            os._exit(1)

    def ngc_matrix(self):
        print("INFO: Using the Rotation matrix for NGC")
        
        axx = ayy = 1 / 2.
        
        axy = - np.sqrt(3) / 2.
        axz = ayz = azx = azy = 0.

        ayx = - axy
        azz = 1



        return [axx, axy, axz, ayx, ayy, ayz, azx, azy, azz]

    def sgc_matrix(self):
        print("INFO: Using the Rotation matrix for SGC")

        axx = ayy = (1 / 2.) + 1. / ( 2. * np.sqrt(2) )
        
        axy = ayx = (1 / 2.) - 1. / ( 2. * np.sqrt(2) )

        axz = - 1. / 2.
        ayz =   1. / 2.
        azx =   1. / 2.
        azy = - 1. / 2.

        azz = 1. / ( np.sqrt(2) )



        return [axx, axy, axz, ayx, ayy, ayz, azx, azy, azz]