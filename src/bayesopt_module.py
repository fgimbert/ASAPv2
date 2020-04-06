from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate
import ase


import time
from pymatgen import MPRester, Composition, Element, Molecule
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
# Import the neccesary tools to generate surfaces
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import sys, getopt
import importlib
import json
import paramiko

import remote_module
importlib.reload(remote_module)


compatibility = MaterialsProjectCompatibility()

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2019, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "November 2019"


domain_dico = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                       'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                       'h': {'name': 'h', 'type': 'continuous', 'domain': (1, 5)},
                       'alpha': {'name': 'alpha', 'type': 'continuous', 'domain': (-180, 180)},
                       'beta': {'name': 'beta', 'type': 'continuous', 'domain': (-180, 180)},
                       'gamma': {'name': 'gamma', 'type': 'continuous', 'domain': (-180, 180)}}


class Bayesopt(object):

    def __init__(self, password, passphrase, MAPI_KEY='4oBTKz0pkFSg9EUQ'):

        self.password = password
        self.passphrase = passphrase

        with open('boinput/cluster.json', mode='r', encoding='utf-8') as fd:
            json_file = json.load(fd)
            self.kwargs = json_file['kwargs']
            self.host = json_file['host']

        if 'passphrase' in list(self.kwargs.keys()):
            self.kwargs['passphrase'] = self.passphrase
        if 'password' in list(self.kwargs.keys()):
            self.kwargs['password'] = self.password

        with open("testfile.txt", "w") as testfile:
            testfile.write('Bonjour, test here\n')

        with open("testoutput.txt", "w") as file:
            file.write('Bonjour, subprocess here\n')
            file.write('\n')

            # slab = ase.io.read('temp/slab', format='vasp')
            file.write('Slab\n')

            with open('boinput/slab', 'r') as slab:
                lines = slab.readlines()
                for line in lines:
                    file.write(line)
            file.write('Molecule\n')
            with open('boinput/molecule', 'r') as slab:
                lines = slab.readlines()
                for line in lines:
                    file.write(line)
            file.write('Domain\n')

            with open('boinput/domain', 'r') as slab:
                lines = slab.readlines()
                Adsorbate = str(line[0])
                file.write(Adsorbate)
                file.write('\n')

                for line in lines[1:]:
                    file.write(line)

            file.write('\n')
            file.flush()

            molecule = ase.io.read('boinput/molecule', format='xyz')
            slab = ase.io.read('boinput/slab', format='vasp')

            file.write('Reading finished\n')

            self.structure = self.add_adsorbate(slab=slab, molecule=molecule, x=0, y=0, h=2.0)

            file.write('Add_adsorbate finished\n')
            ase.io.write('boinput/structure', self.structure, format='vasp')

            # time.sleep(60)

            #self.remote_workdir = 'TestBO/'
            #self.remote = remote_module.Remote(password=self.password, passphrase=self.passphrase, gui=False)
            #self.remote.create_workdir(path=self.remote_workdir)

            # self.remote.put_files(path=self.remote_workdir, path_input=path_input)
            #self.checkls()

            file.write('Finished\n')

    def add_adsorbate(self, slab=None, molecule=None, h=None, x=None, y=None, offset=None):

        with open("testadds.txt", "w") as file:
            file.write('Bonjour, add_adsorbate here\n')
            file.flush()

            from pymatgen.io.ase import AseAtomsAdaptor

            from ase.build import add_adsorbate

            molecule = molecule
            structure = slab.copy()

            file.write('Hello, get_atoms start: \n')
            file.flush()
            # ads_slab = AseAtomsAdaptor.get_atoms(structure)
            cell = structure.get_cell()
            file.write('Hello, get_cell finished\n')

            if h is None:
                h = 2.0

            if x is None and y is None:
                x_cart = 0
                y_cart = 0
            else:
                x_cart = x * cell[0][0] + y * cell[1][0]
                y_cart = x * cell[0][1] + y * cell[1][1]

            # print(h, x, y)

            add_adsorbate(structure, molecule, h, (x_cart, y_cart))
            file.write('Hello, add_adsorbate finished\n')

            return structure

    def checkls(self, path=''):

        """
        Connect remote host and check remote files
        """
        with open("testls.txt", "w") as file:
            with paramiko.SSHClient() as client:
                client.load_system_host_keys()
                client.set_missing_host_key_policy(paramiko.WarningPolicy())

                client.connect(self.host, **self.kwargs)
                stdin, stdout, stderr = client.exec_command("ls {0}".format(path))
                errs = stderr.read().decode('utf-8')

                eerrs = stderr.read().decode('utf-8')
                dir_client = stdout.read().decode('utf-8')
                if eerrs:
                    print(f"ls failed in connecthost:")
                file.write(f"File list: {'  '.join(dir_client.split())}")
                client.close()
