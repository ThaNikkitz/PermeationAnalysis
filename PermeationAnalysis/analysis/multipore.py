import pdb
import prody
import .msd_pf
import .permeation_events
import .axial_loads
import .survival_time

from datetime import datetime

class Multipore:

    def __init__(self, 
                 dcd=None, 
                 pdb=None, 
                 psf=None, 
                 sel='name OH2', 
                 ref='protein and name CA', 
                 upper_check=None, 
                 lower_check=None, 
                 check_radii=None, 
                 delta=3):
        
        self._dcd = dcd
        self._pdb = pdb
        self._psf = psf
        self._sel = sel
        self._ref = ref
        self._upper_check = upper_check
        self._lower_check = lower_check
        self._check_radii = check_radii
        self._delta = delta

    @property
    def dcd(self):
        return self._dcd

    @dcd.setter
    def dcd(self, dcd):
        #TODO: modify so it takes prody trajectory files as input as well
        if isinstance(dcd, list):
            if all([type(element) == str for element in dcd]):
                self._dcd += dcd
            else:
                raise TypeError('The elements of the list must all be strings, i.e. paths to the .dcd files.')
        elif isinstance(dcd, str):
            self._dcd.append(dcd)

        elif isinstance(dcd, (prody.DCDFile, prody.Trajectory)):
            pass

        else:
            raise TypeError('Must be a string pointing to the .dcd file.')
        
    @property
    def pdb(self):
        return self._pdb
    
    @pdb.setter
    def pdb(self, pdb):
        if isinstance(pdb, (prody.AtomGroup, prody.Ensemble, prody.Atomic)):
            self._pdb = pdb
        elif isinstance(pdb, str):
            self._pdb = prody.parsePDB
        else:
            raise TypeError('The input argument is neither a compatible ProDy pdb, nor a path to one such file')

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



