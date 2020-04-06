from ipywidgets import widgets
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen import MPRester
from IPython.display import Image, display, clear_output
from pyscript.localfiles import LocalFiles
from traitlets import traitlets

# from SelectFilesWidget import *
from ipywidgets import Layout
from ase.calculators.vasp import Vasp

import json
import paramiko
from scp import SCPClient
import os

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2019, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "October 2019"


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class configVASP(object):

    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ', VASP_PP_PATH="/home/fgimbert/Projects/potvasp"):
        """Initialize remote class.
           Arguments
           ---------

        """

        self.rester = MPRester(MAPI_KEY)
        self.VASP_PP_PATH = VASP_PP_PATH
        self.accuracy = None
        self.__create_config()

    def __create_config(self):
        self.xc_buttons = widgets.ToggleButtons(options=['LDA', 'PBE'],
                                                value='PBE',
                                                description='xc:',
                                                disabled=False)

        self.setup_buttons = widgets.ToggleButtons(options=['VASP', 'Materials Project'],
                                                   value='Materials Project',
                                                   description='PP choice:',
                                                   tooltips=['corresponds to the table of recommended PAW setups '
                                                             'supplied by the VASP developers.',
                                                             'corresponds to the Materials Project recommended PAW '
                                                             'setups.'],
                                                   disabled=False)

        self.accuracy_buttons = widgets.ToggleButtons(options=['Low', 'Medium', 'High', 'Normal', 'Accurate', 'Single'],
                                                      value='Normal',
                                                      description='PREC:',
                                                      disabled=False)

        self.settings_text = widgets.Textarea(value=' ISIF = 0\n NSW = 0\n NELM = 99\n ISTART = 0',
                                              placeholder='Your own settings',
                                              description='VASP settings:',
                                              disabled=False,
                                              layout={'height': '150px'})

        self.jobfile_text = widgets.Textarea(value='#!/bin/sh\n#$ -cwd\n#$ -N JOBNAME\n#$ -o output\n#$ -V\n#$ -S /bin/bash\n#$ -pe x40 40\nsource /opt/intel/impi/2018.2.199/intel64/bin/mpivars.sh\nexport I_MPI_FALLBACK_DEVICE=0\nexport I_MPI_PIN=0\nexport I_MPI_MPD_RSH=/usr/bin/rsh\nexport I_MPI_HYDRA_BOOTSTRAP=rsh\n\nmpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /usr/local/vasp.5.3.5_intelmpi/bin/vasp ',
                                              placeholder='Job file',
                                              description='',
                                              disabled=False,
                                              layout={'height': '250px', 'width':'100%'})

    def create_input(self, atoms=None, path='input_vasp/', xc='PBE', setup='materialsproject', prec='Accurate'):

        import numpy as np
        # some functions
        def norm(vec):
            v2 = np.dot(vec, vec)
            return np.sqrt(v2)

        def get_kpts(cell, pitch=0.1):
            a1 = cell[0]
            a2 = cell[1]
            a3 = cell[2]
            b1 = 2.0 * np.pi * np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
            b2 = 2.0 * np.pi * np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
            b3 = 2.0 * np.pi * np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
            bb1 = norm(b1)
            bb2 = norm(b2)
            bb3 = norm(b3)
            kpts = (int(bb1 / pitch), int(bb2 / pitch), int(bb3 / pitch))
            return kpts

        from ase import Atoms

        print('Hello')
        if atoms is None:
            a = 4.05  # Gold lattice constant
            b = a / 2
            atoms = Atoms('Li', cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)

        # print(r"C:\Users\f-gimbert\Documents\potvasp")

        # make kpoints with some pitch, which user must specify
        pitch = 2.0 * np.pi / 4.041 / 24
        kpts = get_kpts(atoms.cell, pitch)

        #os.environ['VASP_PP_PATH'] = "/home/fgimbert/Projects/potvasp"
        os.environ['VASP_PP_PATH'] = self.VASP_PP_PATH
        calc = Vasp(xc=xc, setups=setup, kpts=[7, 7, 1])
        calc.initialize(atoms)
        calc.set(prec=prec, ediff=1E-5)
        originpath = os.getcwd()
        os.chdir(path)

        calc.write_incar(atoms)
        calc.write_kpoints()
        calc.write_potcar()

        incar = open('INCAR', 'a', newline='\n')
        incar.write(self.settings_text.value)
        incar.close()

        jobfile = open('submit_file.sh', 'w', newline='\n')
        jobfile.write(self.jobfile_text.value)
        jobfile.close()



        from ase.io import read, write
        write('POSCAR', atoms, format="vasp")

        
        os.chdir(originpath)

    def bulk_input(self, mpid=None, path='input_vasp/', xc='PBE', setup='materialsproject', prec='Accurate'):

        #os.environ['VASP_PP_PATH'] = "/home/fgimbert/Projects/potvasp"
        os.environ['VASP_PP_PATH'] = self.VASP_PP_PATH
        from pymatgen.io.ase import AseAtomsAdaptor

        structure = self.rester.get_structure_by_material_id(mpid)
        atoms = AseAtomsAdaptor.get_atoms(structure)
        calc = Vasp(xc=xc, setups=setup)
        calc.initialize(atoms)
        calc.set(prec=prec, ediff=1E-6)
        originpath = os.getcwd()
        os.chdir(path)

        calc.write_incar(atoms)
        calc.write_kpoints()
        calc.write_potcar()

        from ase.io import read, write
        write('POSCAR', atoms, format="vasp")
        os.chdir(originpath)



