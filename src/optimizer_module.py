from ipywidgets import widgets
from IPython.display import Image, display, clear_output
from traitlets import traitlets
from ipywidgets import Layout

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate

from pymatgen import MPRester, Composition, Element, Molecule
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
# Import the neccesary tools to generate surfaces
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

compatibility = MaterialsProjectCompatibility()

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2019, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "August 2019"


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))
        
        
class Optimizer(object):
    
    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ'):

        self.rester = MPRester(MAPI_KEY)
        
        h = 3.6
        d = 1.10
        slab = fcc111('Al', size=(4, 4, 2), vacuum=10.0)
        molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])
        add_adsorbate(slab, molecule, h, 'ontop')

        self.slab = None
        
        self.calculator = 'EMT'
        self.ninit = 5
        self.__create_runopt_button()
        self.__create_testopt_button()
        self.__create_runcalc_button()
        self.__create_inputvasp_button()
        self.__create_posx_text()
        self.__create_posy_text()
        self.__create_calculator_dropdown()
        self.__create_parameters_select()
        self.__create_calculator_dropdown()
        self.__create_maxbo_slider()
        self.__create_bo_progress()
        self._create_initbo_text()
        self._create_maxbo_text()
        self.__create_initbo_slider()

    def __create_posx_text(self):
        """Build the Text widget for input position."""
        # create widget
        
        self.posx_text = widgets.Text(placeholder='Type something',
                                      value = None,
                                      description='x coordinate',
                                      disabled=False,
                                      style={'description_width': 'initial'})

    def __create_posy_text(self):
        """Build the Text widget for input position."""
        # create widget
        
        self.posy_text = widgets.Text(placeholder='Type something',
                                      value = None,
                                      description='y coordinate',
                                      disabled=False,
                                      style={'description_width': 'initial'})

    def __create_calculator_dropdown(self):
        """Build the dropdown widget to select calculator."""
        
        self.calculator_dropdown = widgets.Dropdown(options=['EMT', 'VASP', 'VASP_inwork'],
                                                    value='EMT')

    def __create_parameters_select(self):
        """Build the SelectMultiple widget to select parameters to be optimized"""
        
        self.parameters_select = widgets.SelectMultiple(options=['x', 'y', 'h'],
                                                        value=['x', 'y', 'h'],
                                                        description='Parameters')

    def _create_maxbo_text(self):
        self.maxbo_text = widgets.HTML(value="Maximum Optimization steps: ")

    def __create_maxbo_slider(self):
        """Build the slider widget to choose Max BO steps."""
        
        self.maxbo_slider = widgets.IntSlider(value=15,
                                              min=5,
                                              max=300,
                                              step=5,
                                              description='')

    def _create_initbo_text(self):
        self.initbo_text = widgets.HTML(value="Number of initial steps:")

    def __create_initbo_slider(self):
        """Build the slider widget to choose Max BO steps."""

        self.initbo_slider = widgets.IntSlider(value=5,
                                              min=5,
                                              max=30,
                                              step=5,
                                              description='')

    def __create_bo_progress(self):
        """Build the progress widget to show  BO steps."""
        maxbo = int(self.maxbo_slider.value)+self.ninit
        self.bo_progress = widgets.IntProgress(value=0,
                                               min=0,
                                               max=maxbo,
                                               step=1,
                                               description='BO Progress: ')
    
    def __create_runcalc_button(self):
        """Build the button widget to run a calculation."""

        self.runcalc_button = LoadedButton(description='Run Calculation',
                                           value=False,
                                           disabled=False,
                                           layout=Layout(width='150px'),
                                           style={'button_color': 'lightblue'})

    def __create_inputvasp_button(self):
        """Build the button widget to create input vasp *in test*."""

        self.inputvasp_button = LoadedButton(description='Create Input',
                                             value=False,
                                             disabled=False,
                                             layout=Layout(width='150px'),
                                             style={'button_color': 'lightblue'})

    def __create_runopt_button(self):
        """Build the button widget to run bayesian optimization."""
        self.runopt_button = LoadedButton(description='Run Optimization',
                                          value=False,
                                          disabled=False,
                                          layout=Layout(width='150px'),
                                          style={'button_color': 'lightblue'})

    def __create_testopt_button(self):
        """Build the button widget to test bayesopt scropt (run only one step)."""
        self.testopt_button = LoadedButton(description='Test Opt. step',
                                          value=False,
                                          disabled=False,
                                          layout=Layout(width='150px'),
                                          style={'button_color': 'lightblue'})

    def calculate(self, repertory='', method='EMT'):
        
        from pymatgen.io.ase import AseAtomsAdaptor
        # pylint: disable=some-message,another-one
        structure = self.runcalc_button.value
        self.slab = structure

        if method == 'EMT':
            self.slab.set_calculator(EMT())
        
            self.e_slab = self.slab.get_potential_energy()
            
        elif method == 'VASP':

            print('VASP not implemented yet. EMT used.')
            self.slab.set_calculator(EMT())
            self.e_slab = self.slab.get_potential_energy()

        elif method == 'VASP_inwork':
            # All steps
            # 1/ Create the input VASP inside work repertory
            # 2/ Send the files on cluster
            # 3/ Run the calculation on cluster
            # 4/ Check the queue every 5mn
            # 5/ If finished, download the ouput files and read energy

            # START
            # 1/ Create the input VASP

            from ase.io import read, write
            write('POSCAR', self.slab, format="vasp")
            self.e_slab = 0

            #self.slab.set_calculator(EMT())
            #self.e_slab = self.slab.get_potential_energy()

    def process_search(self):        
        """Start the Bayesian Optimization Process"""
        print('Bayesian Optimization for {0} adosrbate on {1} surface'.format(1, 2))
