from abc import ABC
from numpy import isin
import prody

import warnings # For checkpoints position check (lower > upper)

class Pore(ABC):
    """
    """

    def __init__(self, 
                 dcd, 
                 pdb, 
                 upper_check,
                 lower_check,
                 sel='name OH2', 
                 ref='protein and name CA'
                 ):
        
        self.dcd = dcd
        self.pdb = pdb
        self.sel = sel
        self.ref = ref
        self.upper_check = upper_check
        self.lower_check = lower_check

    @property
    def dcd(self):
        return self._dcd

    @dcd.setter
    def dcd(self, dcd):
        #TODO: modify so it takes prody trajectory files as input as well
        if isinstance(dcd, list):
            if all([type(element) == str for element in dcd]):
                self._dcd = prody.Trajectory('trajectory')
                for file in dcd:
                    self._dcd.addFile(file)
            else:
                raise TypeError('The elements of the list must all be strings, i.e. paths to the .dcd files.')
            
        elif isinstance(dcd, str):
            self._dcd = prody.DCDFile(dcd)

        elif isinstance(dcd, (prody.DCDFile, prody.Trajectory)):
            self._dcd = dcd

        else:
            raise TypeError('Must be a string pointing to the .dcd file, or a prody DCD or Trajectory.')
        
    @property
    def pdb(self):
        return self._pdb
    
    @pdb.setter
    def pdb(self, pdb):
        if isinstance(pdb, str):
            self._pdb = prody.parsePDB(pdb)
        elif isinstance(pdb, (prody.AtomGroup, prody.Ensemble, prody.Atomic)):
            self._pdb = pdb
        else:
            raise TypeError('The input argument is neither a compatible ProDy pdb, nor a path to one such file')
    
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
    def lower_check(self):
        return self._lower_check
    
    @lower_check.setter
    def lower_check(self, lower_check):
        if isinstance(lower_check, (float, int, str)):
            self._lower_check = lower_check
        else:
            raise TypeError('The upper and lower checkpoints have to be either\
                             a Z position or an atomic selection')

    @property
    def sel(self):
        return self._sel
    
    @sel.setter
    def sel(self, sel):
        if isinstance(sel, str):
            self._sel = self
        else:
            raise TypeError('The selection must be a string')
        
    @property
    def ref(self):
        return self._ref
    
    @ref.setter
    def ref(self, ref):
        if isinstance(ref, str):
            self._ref = ref
        else:
            raise TypeError('The reference selection has to be a string')

    def setup(self):
        self.dcd.link(self.pdb)
        self.dcd.setCoords(self.pdb)
        self.dcd.setAtoms(self.pdb.select(self.ref))
        self.sel = self.pdb.select(self.sel)

        # TODO: Add logic so that only a one-atom-per-residue selection is possible in "sel"
        #       Add logic so as to avoid empty ref and sel selections

        if isinstance(self.upper_check, str):
            self.upper_check = lambda: prody.calcCenter(self.pdb.select(self.upper_check))
        else:
            self.upper_check = lambda: self.upper_check

        if isinstance(self.lower_check, str):
            self.lower_check = lambda: prody.calcCenter(self.pdb.select(self.lower_check))
        else:
            self.lower_check = lambda: self.lower_check

        if isinstance(self.lower_check, (int, float)) and isinstance(self.upper_check, (int, float)):
            if self.lower_check > self.upper_check:
                warnings.warn('The lower checkpoint shouldn\'t be at a higher Z coordinate\
                              than the upper checkpoint. Switching their values')
                dummy = self.lower_check
                self.lower_check = self.upper_check
                self.upper_check = dummy

    def run(self, wrap=True):
        pass
    