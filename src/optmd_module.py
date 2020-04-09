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

from src import remote_module
#importlib.reload(remote_module)

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2019, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "April 2020"



class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class OptMD(object):

    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ', VASP_PP_PATH="/home/fgimbert/Projects/potvasp", ESPRESSO_PSEUDO="/home/fgimbert/Projects/potvasp"):
        """Initialize remote class.
           Arguments
           ---------

        """

        self.rester = MPRester(MAPI_KEY)
        self.VASP_PP_PATH = VASP_PP_PATH
        self.accuracy = None
        self.__create_config()

       

    def __create_config(self):
        
        self.runoptmd_button = LoadedButton(description='Run on-the-fly Opt MD', value=None, disabled=False)

        self.workdir_text = widgets.Text(    
            placeholder='/home/f-gimbert/TESTOptMD',
            value = '/home/f-gimbert/TESTOptMD',
            description='Remote workdir',
            disabled=False,
            style={'description_width': 'initial'}
        )
        
        self.vasppp_text = widgets.Text(    
            placeholder='/home/f-gimbert/potvasp',
            value = '/home/f-gimbert/potvasp',
            description='Vasp PP path',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.ncpus_text = widgets.Text(    
            placeholder='80',
            value = '80',
            description='Number of CPUS',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.vasprun_text = widgets.Text(    
            placeholder='mpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /home/share/vasp.5.4.4/bin/vasp_std',
            value = 'mpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /home/share/vasp.5.4.4/bin/vasp_std',
            description='Vasp PP path',
            disabled=False,
            style={'description_width': 'initial'}
        )

    def create_remote(self):

        self.remote = remote_module.Remote(password=self.password, passphrase=self.passphrase, gui=False)

        with open('../asapdata/cluster.json', mode='r', encoding='utf-8') as fd:
            json_file = json.load(fd)
            self.kwargs = json_file['kwargs']
            self.host = json_file['host']

        if 'passphrase' in list(self.kwargs.keys()):
            self.kwargs['passphrase'] = self.passphrase
        if 'password' in list(self.kwargs.keys()):
            self.kwargs['password'] = self.password

        self.remote.create_workdir(path=self.workdir_text.value, host=self.host, kwargs=self.kwargs)
        print('remote created\n')


    def input_files(self):

        with open('input_optmd/pyRun_script.sh', 'r') as jobfile:
            filedata = jobfile.read()

            # Replace the target string
            
            filedata = filedata.replace('NCPUS', str(int(self.ncpus_text.value)))

            # Write the file out again
            with open('input_optmd/pyRun.sh', 'w', newline='\n') as filejob:
                filejob.write(filedata)

        with open('input_optmd/vasp_setup_script.py', 'r') as jobfile:
            filedata = jobfile.read()

            # Replace the target string
            
            filedata = filedata.replace('VASPPPPATH','"{0}"'.format(self.vasppp_text.value))
            filedata = filedata.replace('ASEVASPCOMMAND', '"{0}"'.format(self.vasprun_text.value))

            # Write the file out again
            with open('input_optmd/vasp_setup.py', 'w', newline='\n') as filejob:
                filejob.write(filedata)

        self.remote.put_files(path=self.workdir_text.value, path_input='input_optmd/', host=self.host, kwargs=self.kwargs)
        print('files moved\n')

    def runscript(self):

        self.remote.execute_command(command='qsub pyRun.sh', path=self.workdir_text.value, host=self.host, kwargs=self.kwargs)

    def runoptmd_action(self, atoms, password, passphrase):

        self.password = password
        self.passphrase = passphrase

        print('in run on the fly aimd')
        self.create_remote()
        self.input_files()
        self.runscript()


        # Create workdir on remote cluster