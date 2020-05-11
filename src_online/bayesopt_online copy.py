from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate
import ase
import os

import time
from pymatgen import MPRester, Composition, Element, Molecule
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
# Import the neccesary tools to generate surfaces
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor

import GPyOpt
from GPyOpt.methods import BayesianOptimization

import sys, getopt
import json
#import paramiko
import numpy as np


from fireworks import LaunchPad

compatibility = MaterialsProjectCompatibility()

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2020, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "April 2020"


domain_dico = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                       'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                       'h': {'name': 'h', 'type': 'continuous', 'domain': (0.7, 4)},
                       'alpha': {'name': 'alpha', 'type': 'continuous', 'domain': (-180, 180)},
                       'beta': {'name': 'beta', 'type': 'continuous', 'domain': (-180, 180)},
                       'gamma': {'name': 'gamma', 'type': 'continuous', 'domain': (-180, 180)}}

domain_dico_discrete = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                                'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                                'h': {'name': 'h', 'type': 'continuous', 'domain': (0.7, 4)},
                                'alpha': {'name': 'alpha', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'beta': {'name': 'beta', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'gamma': {'name': 'gamma', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)}}


class Bayesopt(object):

    def __init__(self, workdir, maxbo=15, initbo=5, x=0.0, y=0.0, MAPI_KEY='4oBTKz0pkFSg9EUQ'):

        self.system = 'Atom'

        self.bo_step = 0
        self.initbo = int(initbo)
        self.maxbo = int(maxbo)
        self.vacuum = 10.0

        self.xads = float(x)
        self.yads = float(y)

        self.emin = None
        self.slabmin = None

        self.workdir = workdir

        self.step_path = ''
        # Construct the domain optimization
        self.domain = []

        #self.localpath = '../asapdata/' + self.workdir
        #os.makedirs(self.localpath, exist_ok=True)

        
        with open("testoutput.txt", "a+") as file:
            file.write('Bonjour, subprocess here\n')
            file.write('\n')

            # slab = ase.io.read('temp/slab', format='vasp')
            file.write('Slab\n')

            with open('input_structure/slab', 'r') as slab:
                lines = slab.readlines()
                for line in lines:
                    file.write(line)
            file.write('Molecule\n')
            with open('input_structure/molecule', 'r') as slab:
                lines = slab.readlines()
                for line in lines:
                    file.write(line)
                    
            file.write('Domain\n')

            with open('input_structure/domain', 'r') as slab:
                lines = slab.readlines()
                Adsorbate = str(line[0])
                file.write(Adsorbate)
                file.write('\n')

                for line in lines[1:]:
                    file.write(line)
                    self.domain.append(domain_dico_discrete[str(line.strip('\n'))])

            file.write('\n')
            file.flush()

            self.molecule = ase.io.read('input_structure/molecule', format='xyz')
            self.slab = ase.io.read('input_structure/slab', format='vasp')

            file.write('Reading finished\n')

            self.structure = self.add_adsorbate_ase(slab=self.slab, molecule=self.molecule,
                                                    x=self.xads, y=self.yads, h=2.0)

            file.write('Add_adsorbate finished\n')

            ase.io.write('structure', self.structure, format='vasp')
            ase.io.write('/input_vasp/POSCAR', self.structure, format='vasp')

            # time.sleep(60)

            #self.remote_workdir = 'ASAP/BO/' + self.workdir
            #self.remote_input = self.remote_workdir + "/input"
            #self.remote = remote_module.Remote(password=self.password, passphrase=self.passphrase, gui=False)
            #file.write('remote created\n')
            #self.remote.create_workdir(path=self.remote_input, host=self.host, kwargs=self.kwargs)
            #file.write('remote created\n')

            #self.remote.put_files(path=self.remote_input, path_input=self.localpath+'/vasp/', host=self.host, kwargs=self.kwargs)
            #file.write('files moved\n')

            # self.remote.put_files(path=self.remote_workdir, path_input=path_input)

            #file.write('checkls here\n')
            #file.write('\n')
            #self.checkls()
            #file.write('checkls finished\n')

            with open('input_vasp/submit_file.sh', 'r') as jobfile:
                self.submitfile = jobfile.read()

            file.write('Start BO')
            #self.runopt()

            outf = open("output.txt", "a+")
            outf.write('Finished')
            outf.close()
            file.write('Finished\n')

            

    def add_adsorbate_ase(self, slab=None, molecule=None, h=None, x=None, y=None, offset=None):

        with open("../asapdata/{0}/testadds.txt".format(str(self.workdir)), "a+") as file:
            file.write('Hello, add_adsorbate here\n')
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
                x_cart = self.xads
                y_cart = self.yads
            else:
                x_cart = x * cell[0][0] + y * cell[1][0]
                y_cart = x * cell[0][1] + y * cell[1][1]

            # print(h, x, y)

            add_adsorbate(structure, molecule, h, (x_cart, y_cart))
            structure.center(axis=2, vacuum=self.vacuum)
            file.write('Hello, add_adsorbate finished\n')

            return structure

    def add_adsorbate_ase_mol(self, slab=None, molecule=None, h=None, x=None, y=None, alpha=None, beta=None, gamma=None,
                              offset=None):

        with open("../asapdata/{0}/testadds.txt".format(str(self.workdir)), "a+") as file:
            file.write('Hello, add_adsorbate here\n')
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

            if alpha is None:
                molecule.rotate(0, 'x')
            else:
                # print(type(alpha))
                alpha = float(alpha)
                # print(type(alpha))
                molecule.rotate(alpha, 'x')
                # print('rotation done')

            if beta is None:
                molecule.rotate(0, 'y')
            else:
                beta = float(beta)
                molecule.rotate(beta, 'y')

            if gamma is None:
                molecule.rotate(0, 'z')
            else:
                gamma = float(gamma)
                molecule.rotate(gamma, 'z')

            if h is None:
                h = 3.0

            if x is None and y is None:
                x_cart = self.xads
                y_cart = self.yads
            else:
                x_cart = x * cell[0][0] + y * cell[1][0]
                y_cart = x * cell[0][1] + y * cell[1][1]

            # print(h, x, y)

            add_adsorbate(structure, molecule, h, (x_cart, y_cart))
            structure.center(axis=2, vacuum=self.vacuum)
            file.write('Hello, add_adsorbate finished\n')

            return structure

    def checkls(self, path='ASAP/BO/'):

        """
        Connect remote host and check remote files
        """

        with open("../asapdata/{0}/testls.txt".format(str(self.workdir)), "a+") as file:
            with paramiko.SSHClient() as client:
                client.load_system_host_keys()
                client.set_missing_host_key_policy(paramiko.WarningPolicy())

                client.connect(self.host, **self.kwargs)
                stdin, stdout, stderr = client.exec_command("ls {0}".format(self.remote_workdir))
                errs = stderr.read().decode('utf-8')

                eerrs = stderr.read().decode('utf-8')
                dir_client = stdout.read().decode('utf-8')
                if eerrs:
                    print(f"ls failed in connecthost:")

                for line in dir_client.split():
                    file.write(line)
                    file.write('\n')

                client.close()

    def runopt(self):

        # For 6d optimizations, it would be better to do several optimizations than a big one. Indeed, a big
        # optimizations can be stucked in a low minima trying to find the best angles far away from the surface.

        # THe idea is to do:
        # 1/ free optimizations 10 - 15 steps but with large angle value steps (90 for example)
        # 2/ take the best result and optimize only the x,y,h parameters.
        # 3/ take the best and optimize the angle with better accuracy (free or 10 degrees step)
        # 4/ free optimization (probably useless, the BO will explore again the angles values)

        # Another possibility is to modify the acquisition function to improve the exploitation and not the
        # exploration. by default, jitter = 0.01. Increase improve exploration, decrease exploitation.
        # Try several parameters
        # Calculate the next step with 2 different BO, one with normal parameter, one with more exploitation
        # Each time select randomly the next step between the 2 values

        myBopt = BayesianOptimization(f=self.bo_calculate, domain=self.domain,
                                      initial_design_numdata=self.initbo,
                                      acquisition_type='EI',
                                      exact_feval=True)

        myBopt.run_optimization(max_iter=self.maxbo)
        myBopt.plot_convergence(filename="../asapdata/{0}/bo_convergence".format(str(self.workdir)))
        if len(self.domain) <=2:
            myBopt.plot_acquisition(filename="../asapdata/{0}/bo_acquisition".format(str(self.workdir)))

        # print(x[0])
        x_input = None
        y_input = None
        h_input = None
        alpha_input = None
        beta_input = None
        gamma_input = None
        # print(self.domain)

        for i, parameter in enumerate(self.domain):
            if parameter['name'] == 'x':
                x_input = myBopt.x_opt[i]
            if parameter['name'] == 'y':
                y_input = myBopt.x_opt[i]
            if parameter['name'] == 'h':
                h_input = myBopt.x_opt[i]
            if parameter['name'] == 'alpha':
                alpha_input = myBopt.x_opt[i]
            if parameter['name'] == 'beta':
                beta_input = myBopt.x_opt[i]
            if parameter['name'] == 'gamma':
                gamma_input = myBopt.x_opt[i]

    def bo_calculate(self, x):

        with open("../asapdata/{0}/testbo.txt".format(str(self.workdir)), "a+") as file:

            file.write('Read input files\n')

            self.molecule = ase.io.read(self.localpath+'/molecule', format='xyz')
            self.slab = ase.io.read(self.localpath+'/slab', format='vasp')

            self.step_path = self.localpath + '/bo_' + str(self.bo_step)
            file.write(self.step_path)
            file.write(' \n')
            file.flush()

            if not os.path.exists(self.step_path):
                os.makedirs(self.step_path)
            file.write('directory created files\n')
            file.flush()

            x_input = None
            y_input = None
            h_input = None
            alpha_input = None
            beta_input = None
            gamma_input = None

            for i, parameter in enumerate(self.domain):
                if parameter['name'] == 'x':
                    x_input = x[0][i]
                if parameter['name'] == 'y':
                    y_input = x[0][i]
                if parameter['name'] == 'h':
                    h_input = x[0][i]
                if parameter['name'] == 'alpha':
                    alpha_input = x[0][i]
                if parameter['name'] == 'beta':
                    beta_input = x[0][i]
                if parameter['name'] == 'gamma':
                    gamma_input = x[0][i]

                    # structure = self.builder.add_adsorbate_ase(h=x[0][0])
            if self.system == 'Atom':
                self.structure = self.add_adsorbate_ase(slab=self.slab, molecule=self.molecule, x=x_input, y=y_input,
                                                        h=h_input)

            elif self.system == 'Molecule':
                # Molecule not implemented yet !!
                self.structure = self.add_adsorbate_ase_mol(slab=self.slab, molecule=self.molecule, x=x_input,
                                                            y=y_input, h=h_input, alpha=alpha_input, beta=beta_input,
                                                            gamma=gamma_input)

            # Run the VASP calculation for self.bo_step !
            y_next = self.calculate(x)
            print(x, y_next)

            self.bo_step += 1

            return y_next

    def calculate(self, x):

        input_vector = ''

        for i, parameter in enumerate(self.domain):
            if parameter['name'] == 'x':
                x_input = x[0][i]
                input_vector += ' x=' + str(x_input)
            if parameter['name'] == 'y':
                y_input = x[0][i]
                input_vector += ' y=' + str(y_input)
            if parameter['name'] == 'h':
                h_input = x[0][i]
                input_vector += ' h=' + str(h_input)
            if parameter['name'] == 'alpha':
                alpha_input = x[0][i]
                input_vector += ' a=' + str(alpha_input)
            if parameter['name'] == 'beta':
                beta_input = x[0][i]
                input_vector += ' b=' + str(beta_input)
            if parameter['name'] == 'gamma':
                gamma_input = x[0][i]
                input_vector += ' g=' + str(gamma_input)

        with open("log.txt", "a+") as file:

            outf = open("output.txt", "a+")
            file.write('BO Step {0} for structure {1} !\n'.format(self.bo_step, self.workdir))
            file.write('Input vector: {0}\n'.format(x[0]))
            file.flush()

            """
            with open(self.search_path + '/inputs_bo.txt', 'a+') as file:
                file.write(str(x[0]) + '\n')

            with open(self.search_path + '/inputs_bofull.txt', 'a+') as file:
                file.write(str(x[0]) + '\n')

            with open(self.search_path + '/input_vasp/POSCAR_model.txt', 'r') as file:
                poscar = file.readlines()
            """

            # Replace the target string

            with open('input_vasp/submit_file.sh', 'r') as jobfile:
                filedata = jobfile.read()

            # Replace the target string
            jobname = self.workdir + '/' + str(self.bo_step)
            file.write(jobname+'\n')
            file.flush()
            jobname = jobname.replace('/', '_')
            filedata = filedata.replace('JOBNAME', jobname)

            # Write the file out again
            with open(self.step_path + '/submit_file.sh', 'w', newline='\n') as filejob:
                filejob.write(filedata)

            file.write('POSCAR creation\n')
            ase.io.write(self.step_path + '/POSCAR', self.structure, format='vasp')
            file.write('POSCAR created inside bo_step {0} directory\n'.format(self.bo_step))

            import shutil
            #Already in workdir on ATLAS
            shutil.copy2('input_vasp/KPOINTS'.format(str(self.workdir)), self.step_path)  # complete target filename given
            shutil.copy2("input_vasp/INCAR".format(str(self.workdir)), self.step_path)  # complete target filename given
            shutil.copy2("input_vasp/POTCAR".format(str(self.workdir)), self.step_path)  # complete target filename given
            shutil.copy2("input_vaspp/job_submit".format(str(self.workdir)), self.step_path)  # complete target filename given

            #self.remote_step = self.remote_workdir + '/bo_' + str(self.bo_step)

            #self.remote.create_workdir(path=self.remote_step, host=self.host, kwargs=self.kwargs)
            file.write('remote created\n')

            file.write('{0}\n'.format(self.remote_step))
            file.write('{0}\n'.format(self.step_path))

            file.flush()

            #self.remote.put_files(path=self.remote_step, path_input=self.step_path+'/', host=self.host, kwargs=self.kwargs)

            test_vasp = True

            if test_vasp:
                # Execute VASP calculations inside test_i on enemat by ssh

                ## CONNECTION TO ATLAS NEED TO CHANGE THIS TO RUN FW JOB !!!


                with paramiko.SSHClient() as client:
                    client.load_system_host_keys()
                    client.set_missing_host_key_policy(paramiko.WarningPolicy)

                    # Establish SSH connection
                    client.connect(self.host, **self.kwargs)
                    # Submit our Grid Engine job by running a remote 'qsub' command over SSH

                    if test_vasp:
                        file.write('connected to {0}\n'.format(self.remote_step))
                        file.flush()

                        ## CONNECTION TO ATLAS NEED TO CHANGE THIS TO RUN FW JOB !!!
                        stdin, stdout, stderr = client.exec_command("chmod +x {0}/job_submit; cd {0}/; ./job_submit".
                                                                    format(self.remote_step))
                        eerrs = stderr.read().decode('utf-8')
                        if eerrs:
                            file.write(f"chmod failed in connecthost:")
                            file.flush()
                            file.write(eerrs)
                            file.flush()
                        file.write('job submitted to {0}\n'.format(self.bo_step))
                        file.flush()

                        id_job = stdout.readlines()[0].split()[2]
                        file.write('idjob: {0}\n'.format(id_job))
                        file.write("qstat -r | grep -cw '{0}'\n".format(id_job))
                        file.flush()
                        time.sleep(10)
                        stdin, stdout, stderr = client.exec_command("qstat -r | grep -cw '{0}'".format(id_job))
                        result = stdout.readlines()
                    else:
                        result = [0]

                while int(result[0]) != 0:
                    with paramiko.SSHClient() as client:
                        client.load_system_host_keys()
                        client.set_missing_host_key_policy(paramiko.WarningPolicy)
                        # Establish SSH connection
                        client.connect(self.host, **self.kwargs)
                        # stdin,stdout,stderr = client.
                        # exec_command("qstat -r | grep -cw {0}".format(prj_info['prefixes'][0]))
                        stdin, stdout, stderr = client.exec_command("qstat -r | grep -cw '{0}'".format(id_job))
                        result = stdout.readlines()
                    time.sleep(300)

                file.write('VASP calculation finished.\n')
                file.flush()

                # Download output vasp calculation repertory self.step_path inside local workdir
                vaspfile = self.remote_step + '/vasprun.xml'
                #try:
                 #   self.remote.get_files(remote_path=vaspfile, local_path=self.step_path, host=self.host, kwargs=self.kwargs)
                #except:
                    # file.write('Download of Vasp files error\n')
                 #   file.flush()

                file.write('Download of Vasp files done. Vasprun Read.\n')
                # Read Energy from vasprun.xml !!!
                file.write('before vasprun \n')
                file.flush()
                if os.path.exists('{0}/vasprun.xml'.format(self.step_path)):
                    file.write('{0}/vasprun.xml exists !\n'.format(self.step_path))

                from pymatgen.io.vasp import Vasprun
                try:
                    vasprun = Vasprun('{0}/vasprun.xml'.format(self.step_path))
                except:
                    file.write('Error while reading vasprun ! Output energy above hull set to 100 \n')
                    file.write('Probably an error during scf. Check input file for BO step {0}\n'.format(self.bo_step))
                    file.write(
                        '{0} {1}/{2} ({4}) {3} -error {5}\n'.format(self.workdir, self.bo_step, self.maxbo, 100,
                                                                    self.initbo, input_vector))
                    outf.write('{0} {1}/{2} ({4}) {3} -error {5}\n'.format(self.workdir, self.bo_step, self.maxbo, 100,
                                                                           self.initbo, input_vector))
                    file.flush()
                    outf.flush()
                    outf.close()
                    return 100

                file.write('after vasprun\n')
                file.flush()
                if vasprun.converged_electronic:
                    file.write('Electronic step convergence has been reached in the final ionic step\n')
                else:
                    file.write('Electronic step convergence has not been reached in the final ionic step\n')

                energy = vasprun.final_energy
                file.write('{0} {1}/{2} ({4}) {3} {5}\n'.format(self.workdir, self.bo_step, self.maxbo, energy,
                                                                self.initbo, input_vector))
                outf.write(
                    '{0} {1}/{2} ({4}) {3} {5}\n'.format(self.workdir, self.bo_step, self.maxbo, energy, self.initbo,
                                                         input_vector))
                file.flush()
                outf.flush()
                outf.close()
                return energy
            else:
                file.write(
                    '{0} {1}/{2} ({4}) {3} {5}\n'.format(self.workdir, self.bo_step, self.maxbo, self.bo_step,
                                                         self.initbo, input_vector))
                outf.write(
                    '{0} {1}/{2} ({4}) {3} {5}\n'.format(self.workdir, self.bo_step, self.maxbo, self.bo_step,
                                                         self.initbo, input_vector))
                file.write('Output energy: {0}\n'.format(self.bo_step))
                file.flush()
                outf.flush()
                outf.close()
                return self.bo_step


def main(argv):

    #lpad = LaunchPad.auto_load()

    opts, args = getopt.getopt(argv, ":", ["workdir=", "maxbo=", "initbo=", "x=", "y="])

    for opt, arg in opts:
        if opt == "--workdir":
            workdir = arg
        elif opt == "--maxbo":
            maxbo = arg
        elif opt == "--initbo":
            initbo = arg
        elif opt == "--x":
            x = arg
        elif opt == "--y":
            y = arg

    bayesopt = Bayesopt(workdir, maxbo, initbo, x, y)


if __name__ == "__main__":

    
    main(sys.argv[1:])