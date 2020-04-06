from ipywidgets import widgets
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen import MPRester
from IPython.display import Image, display, clear_output
from pyscript.localfiles import LocalFiles
from traitlets import traitlets

# from SelectFilesWidget import *
from ipywidgets import Layout
from ase.calculators.espresso import Espresso

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
__date__ = "March 2020"


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class configQE(object):

    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ',ESPRESSO_PSEUDO="/home/fgimbert/Projects/potvasp"):
        """Initialize remote class.
           Arguments
           ---------

        """

        self.rester = MPRester(MAPI_KEY)
        self.ESPRESSO_PSEUDO = ESPRESSO_PSEUDO
        self.accuracy = None
        self.__create_config()

    def __create_config(self):
        self.xc_buttons = widgets.ToggleButtons(options=['LDA', 'PBE'],
                                                value='PBE',
                                                description='xc:',
                                                disabled=False)

        self.jobfile_text = widgets.Textarea(value='#!/bin/sh\n#$ -cwd\n#$ -N JOBNAME\n#$ -o output\n#$ -V\n#$ -S /bin/bash\n#$ -pe x40 40\nsource /opt/intel/impi/2018.2.199/intel64/bin/mpivars.sh\nexport I_MPI_FALLBACK_DEVICE=0\nexport I_MPI_PIN=0\nexport I_MPI_MPD_RSH=/usr/bin/rsh\nexport I_MPI_HYDRA_BOOTSTRAP=rsh\n\nmpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /usr/local/vasp.5.3.5_intelmpi/bin/vasp ',
                                              placeholder='Job file',
                                              description='',
                                              disabled=False,
                                              layout={'height': '250px', 'width':'100%'})

    def create_input(self, atoms=None, path='input_qe/', xc='PBE', setup='materialsproject', prec='Accurate'):

        def norm(vec):
            v2 = np.dot(vec, vec)
            return np.sqrt(v2)

        def create_pseudopot(self):

            path = self.ESPRESSO_PSEUDO
            list_files = [x for x in os.walk(path)][0][2]





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

        os.environ['ESPRESSO_PSEUDO'] = self.ESPRESSO_PSEUDO

        # make kpoints with some pitch, which user must specify
        pitch = 2.0 * np.pi / 4.041 / 24
        kpts = get_kpts(atoms.cell, pitch)

        calc = Espresso(tstress=True, tprnfor=True, kpts=(3, 3, 3))
        calc.write_input()

