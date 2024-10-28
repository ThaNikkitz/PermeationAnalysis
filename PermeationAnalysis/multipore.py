import pdb
from . import Pore
import prody
import .msd_pf
import .permeation_events
import .axial_loads
import .survival_time

import warnings # For checkpoints position check (lower > upper)

from datetime import datetime

class Multipore(Pore):

    def __init__(self, 
                 dcd=None, 
                 pdb=None, 
                 sel='name OH2', 
                 ref='protein and name CA',
                 upper_check=None,
                 lower_check=None, 
                 radius=5, 
                 delta=3
                 ):
        super().__init__(dcd, pdb, sel, ref)
        # TODO: Implement a way to analyze multiple pores at once
        self.upper_check = upper_check
        self.lower_check = lower_check
        self.check_radii = radius
        self.delta = delta

    @property
    def upper_check(self):
        return self._upper_check
    
    @upper_check.setter
    def upper_check(self, upper_check):
        if isinstance(upper_check, (float, int, str)):
            self._upper_check = upper_check
        else:
            raise TypeError('The upper and lower checkpoints have to be either\
                             a Z position or an atomic selection')
        

    @property
    def upper_check(self):
        return self._upper_check
    
    @upper_check.setter
    def upper_check(self, check):
        if self._lower_check is not None and isinstance(check, (float, int)):
            if check < self._lower_check:
                raise ValueError('The lower checkpoint can\'t be at a higher z than the upper one.')
            else:
                self._upper_check = check
        
        elif not isinstance(check, (float, int)):
            raise TypeError('The z position of the upper checkpoint must be a number')
        
        else:
            self._upper_check = check
        
        
    @property
    def lower_check(self):
        return self._lower_check
    
    @upper_check.setter
    def lower_check(self, check):
        if self._upper_check is not None and isinstance(check, (float, int)):
            if check < self._upper_check:
                raise ValueError('The lower checkpoint can\'t be at a higher z than the upper one.')
            else:
                self._lower_check = check
        
        elif not isinstance(check, (float, int)):
            raise TypeError('The z position of the upper checkpoint must be a number')
        
        else:
            self._lower_check = check

    @property
    def check_radii(self):
        return self._check_radii
    
    @check_radii.setter
    def check_radii(self, radii):
        if isinstance(radii, (float, int)):
            if radii > 0:
                self._check_radii = radii
                
            else:
                raise ValueError('Radius of the checkpoints must be a positive number')

        else:
            raise TypeError('Radius of the checkpoints must be a number')            



    def pf(self, bin_number, bin_size=None, ):
        pass

    def permeation_events(self, ):
        pass

    def axial_loads(self, ):
        pass

    def survival_time(self, threshold=(5, 100, 500)):
        pass

    def run_analysis(self, ):
        pass



