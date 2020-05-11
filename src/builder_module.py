from ipywidgets import widgets
from IPython.display import Image, display, clear_output
from traitlets import traitlets
from ipywidgets import Layout
import numpy as np
from ase import Atoms

from ipyfilechooser import FileChooser

from pymatgen import MPRester, Composition, Element, Molecule
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
# Import the neccesary tools to generate surfaces
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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
        
        
class Builder(object):
    
    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ'):

        self.rester = MPRester(MAPI_KEY)
        
        self.__create_miller_text()
        self.__create_vacuum_text()
        self.__create_slabsize_text()
        self.__create_slabcenter_checkbox()
        self.__create_slab_button()
        self.__create_supercell_text()
        self.__create_orientation_button()
        self.__create_slabs_dropdown()
        self.__create_senergy_button()
        self.__create_all_energies_button()
        self.__create_sites_dropdown()
        self.__create_height_slider()
        self.__create_alpha_slider()
        self.__create_beta_slider()
        self.__create_gamma_slider()
        self.__create_orientation_select()
        self.__create_adatom_radio()
        self.__create_adatom_text()
        self.__create_molad_button()
        self.__create_importmol_button()
        self.__create_uploadmol_button()
        self.__create_molad_dropdown()
        
        self.__create_adsorbate_button()

        self.adatom_radio.observe(self.update_adatom, 'value')

        self.ads_slab = None
        self.entry_id = None
        self.ads_sites = None

    def __create_miller_text(self):
        """Build the Text widget for Miller index  input"""
        self.miller_text = widgets.Text(    
            placeholder='h k l',
            value='0 0 1',
            description='Miller index',
            disabled=False,
            layout=Layout(width='150px'),
            style={'description_width': 'initial'}
        )    
        
    def __create_vacuum_text(self):
        """Build the Text widget for min_vacuum_size  input"""
        self.vacuum_text = widgets.Text(    
            placeholder='10',
            value='10',
            description='min_vacuum_size (in Angstroms)',
            disabled=False,
            style={'description_width': 'initial'}
        )

    def __create_slabsize_text(self):
        """Build the Text widget for min_vacuum_size  input"""
        self.slabsize_text = widgets.Text(    
            placeholder='10',
            value = '10',
            description='min_slab_size (in Angstroms)',
            disabled=False,
            style={'description_width': 'initial'}
        )
        
    def __create_slabcenter_checkbox(self):
        """Build the Checkbox widget for Center Slab."""
        self.slabcenter_checkbox = widgets.Checkbox(value=True,
                                                    description='Center Slab',
                                                    style={'description_width': 'initial'})

    def __create_supercell_text(self):
        """Build the Text widget for min_vacuum_size  input"""
        self.supercellx_text = widgets.Text(    
            placeholder='2',
            value = '2',
            description = 'x',
            disabled=False,
            layout=Layout(width='50px'),
            style={'description_width': 'initial'}
        )
        
        self.supercelly_text = widgets.Text(    
            placeholder='2',
            value='2',
            description='y',
            disabled=False,
            layout=Layout(width='50px'),
            style={'description_width': 'initial'}
        )
        
        self.supercellz_text = widgets.Text(    
            placeholder='1',
            value='1',
            description='z',
            disabled=False,
            layout=Layout(width='50px'),
            style={'description_width': 'initial'}
        )

    def __create_orientation_select(self):
        """Build the Select widget for search results ."""

        # create widget
        self.orientation_select = widgets.Select(options=[],
                                                 value=None,
                                                 rows=3,
                                                 description='')

    def __create_slabs_dropdown(self):
        """Build the Text widget for Miller index  input"""
        self.slabs_dropdown = widgets.Dropdown(options=['None'],
                                               value=None,
                                               description='All orientation slabs',
                                               style={'description_width': 'initial'}
        )

    def __create_orientation_button(self):
        """Build the button widget to create slab."""

        # create widget
        self.orientation_button = LoadedButton(description='Check Orientation', value=False, disabled=False)

    def __create_slab_button(self):
        """Build the button widget to create slab."""
        
        # create widget
        self.slab_button = LoadedButton(description='Generate Slab', value=False, disabled=True, button_style='warning')

    def __create_senergy_button(self):
        """Build the button widget to calculate surface energy."""

        # create widget
        self.sernergy_button = LoadedButton(description='Surface Energy', value=False, disabled=True)

    def __create_all_energies_button(self):
        """Build the button widget to calculate surface energy."""

        # create widget
        self.all_energies_button = LoadedButton(description='All Surface Energies', value=False, disabled=True)

    def __create_adatom_radio(self):
        """Build the Radio widget for display outputs."""
        self.adatom_radio = widgets.RadioButtons(options=['Atom', 'Molecule', 'Import Molecule'],
                                                 value=None,
                                                 layout=Layout(width='150px'),
                                                 disabled=False)
        
    def __create_adatom_text(self):
        """Build the Text widget for adatom Text input"""
        self.adatom_text = widgets.Text(disabled=False,
                                        value=None,
                                        placeholder='H, H2O, C60, ...',
                                        layout=Layout(width='150px'),
                                        style={'description_width': 'initial'})

    def __create_importmol_button(self):
        """Build the Upload CIF button widget."""
        self.importmol = FileChooser()
        # update function
        #self.uploadcif_button.on_click(self.uploadcif_clicked)
    
    def __create_uploadmol_button(self):
        """Build the Upload CIF button widget."""

        self.uploadmol_button = LoadedButton(description='Upload Molecule', value=None, disabled=True, button_style='warning')
        
        # update function

    def __create_molad_button(self):
        """Build the button widget to check Molecule."""
        # create widget
        self.molad_button = LoadedButton(description='Check Molecule', 
                                         value=False,
                                         disabled=True,
                                         layout=Layout(width='150px'),
                                         style={'button_color':'lightblue'})

    def __create_molad_dropdown(self):
        """Build the dropdown widget with results from molad_button."""
        
        # create widget
        self.molad_dropdown = widgets.Dropdown(options=[None],
                                               disabled=True,
                                               style={'description_width': 'initial'})

    def update_adatom(self, *args):
        
        if self.adatom_radio.value == 'Atom':
            self.molad_button.disabled = True
            self.uploadmol_button.disabled = True
            self.adatom_text.value = 'H'
            self.adatom_text.placeholder = 'H, O, Au, ...'
            self.alpha_slider.disabled = True
            self.beta_slider.disabled = True
            self.gamma_slider.disabled = True
            
        elif self.adatom_radio.value == 'Molecule':
            self.molad_button.disabled = False
            self.uploadmol_button.disabled = True
            self.adatom_text.value = 'H2O'
            self.adatom_text.placeholder = 'H2O, C60, ...'
            self.alpha_slider.disabled = False
            self.beta_slider.disabled = False
            self.gamma_slider.disabled = False
        elif self.adatom_radio.value == 'Import Molecule':
            self.molad_button.disabled = True
            self.uploadmol_button.disabled = False
            
        else:
            self.molad_button.disabled = True
            self.uploadmol_button.disabled = True
            self.adatom_text.value = None
            self.adatom_text.placeholder = 'H, H2O, C60, ...'

    def __create_adsorbate_button(self):
        """Build the button widget to create slab."""
        
        # create widget
        self.adsorbate_button = LoadedButton(description='Add Adsorbate', value=False, disabled=True)
        
        # update function
        # self.adsorbate_button.on_click(self.add_adsorbate)

    def __create_sites_dropdown(self):
        """Build the button widget to create slab."""
        
        # create widget
        self.sites_dropdown = widgets.Dropdown(options=[None],
                                               description='Relevant sites adsorption',
                                               style={'description_width': 'initial'})

    def __create_height_slider(self):
        """Build the Height slider widget to control the height of adsorbate interactively."""
        self.height_slider = widgets.FloatSlider(value=2.0,
                                                 min=0.5,
                                                 max=self.vacuum_text.value,
                                                 step=0.1,
                                                 description='h:',
                                                 layout=Layout(width='200px'),
                                                 style={'description_width': 'initial'},
                                                 disabled=False)

    def __create_alpha_slider(self):
        """Build the Height slider widget to control the height of adsorbate interactively."""
        self.alpha_slider = widgets.FloatSlider(value=0,
                                                min=-180,
                                                max=180,
                                                step=5,
                                                description='alpha:',
                                                layout=Layout(width='200px'),
                                                style={'description_width': 'initial'},
                                                disabled=False)
    
    def __create_beta_slider(self):
        """Build the Height slider widget to control the angle beta of adsorbate interactively."""
        self.beta_slider = widgets.FloatSlider(value=0,
                                               min=-180,
                                               max=180,
                                               step=5,
                                               description='beta:',
                                               layout=Layout(width='200px'),
                                               style={'description_width': 'initial'},
                                               disabled=False)

    def __create_gamma_slider(self):
        """Build the Height slider widget to control the angle gamma of adsorbate interactively."""
        self.gamma_slider = widgets.FloatSlider(value=0,
                                                min=-180,
                                                max=180,
                                                step=5,
                                                description='gamma:',
                                                layout=Layout(width='200px'),
                                                style={'description_width': 'initial'},
                                                disabled=False)

    @staticmethod
    def is_structure_invertible(structure):
        '''
        This function figures out whether or not an `pymatgen.Structure` object has
        symmetricity.  In this function, the affine matrix is a rotation matrix
        that is multiplied with the XYZ positions of the crystal. If the z,z
        component of that is negative, it means symmetry operation exist, it could
        be a mirror operation, or one that involves multiple rotations/etc.
        Regardless, it means that the top becomes the bottom and vice-versa, and
        the structure is the symmetric.  i.e. structure_XYZ = structure_XYZ*M.
        Arg:
            structure   A `pymatgen.Structure` object.
        Returns
            A boolean indicating whether or not your `ase.Atoms` object is
            symmetric in z-direction (i.e. symmetric with respect to x-y plane).
        '''
        # If any of the operations involve a transformation in the z-direction,
        # then the structure is invertible.
        sga = SpacegroupAnalyzer(structure, symprec=0.1)
        for operation in sga.get_symmetry_operations():
            xform_matrix = operation.affine_matrix
            z_xform = xform_matrix[2, 2]
            if z_xform == -1:
                return True

        return False

    @staticmethod
    def flip_atoms(structure):
        '''
        Flips an atoms object upside down. Normally used to flip slabs.
        Arg:
            atoms   `ase.Atoms` object
        Returns:
            atoms   The same `ase.Atoms` object that was fed as an argument,
                    but flipped upside down.
        '''

        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms = atoms.copy()

        # This is black magic wizardry to me. Good look figuring it out.
        atoms.wrap()
        atoms.rotate(180, 'x', rotate_cell=True, center='COM')
        if atoms.cell[2][2] < 0.:
            atoms.cell[2] = -atoms.cell[2]
        if np.cross(atoms.cell[0], atoms.cell[1])[2] < 0.0:
            atoms.cell[1] = -atoms.cell[1]
        atoms.wrap()

        structure = AseAtomsAdaptor.get_structure(atoms)
        return structure

    def check_orientation(self, entry_id):

        self.entry_id = entry_id

        if len(self.miller_text.value.split()) > 1:

            miller_index = (int(self.miller_text.value.split()[0]),
                            int(self.miller_text.value.split()[1]),
                            int(self.miller_text.value.split()[2]))

            structure = self.rester.get_structure_by_material_id(entry_id)
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            conventional_structure = SpacegroupAnalyzer(structure). \
                get_conventional_standard_structure(international_monoclinic=True)

            # slabgen = generate_all_slabs(conventional_structure,
            # max_index=1, min_slab_size=int(self.slabsize_text.value),
            # min_vacuum_size=int(self.vacuum_text.value), lll_reduce =True, center_slab = self.slabcenter_checkbox.value)

            # slab_dict = {slab.miller_index:slab for slab in slabgen}

            # self.all_slabs = [slab for slab in slabgen if slab.miller_index==miller_index]

            slabgen = SlabGenerator(conventional_structure, miller_index,
                                    int(self.slabsize_text.value),
                                    int(self.vacuum_text.value), lll_reduce=True,
                                    center_slab=self.slabcenter_checkbox.value, max_normal_search=5)

            self.all_slabs = slabgen.get_slabs()

            list_orientations = []
            self.full_slabs = []
            i = 0

            for slab in self.all_slabs:
                orientation = str(i) + '- ' + str(slab.miller_index[0]) + ' ' + str(slab.miller_index[1]) + ' ' \
                              + str(slab.miller_index[2])

                list_orientations.append(orientation)
                self.full_slabs.append(slab)
                i += 1

                if not self.is_structure_invertible(slab):
                    orientation = str(i) + '- ' + str(slab.miller_index[0]) + ' ' + str(slab.miller_index[1]) + ' ' \
                                   + str(slab.miller_index[2]) + ' flip'
                    list_orientations.append(orientation)

                    slab = self.flip_atoms(slab)
                    self.full_slabs.append(slab)
                    i += 1

            return self.full_slabs, len(self.full_slabs), list_orientations, self.all_slabs
        else:  # Max Miller index (for example, 3 for (1, 1, 1))

            structure = self.rester.get_structure_by_material_id(entry_id)
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            conventional_structure = SpacegroupAnalyzer(structure). \
                get_conventional_standard_structure(international_monoclinic=True)

            self.all_slabs = generate_all_slabs(conventional_structure, int(self.miller_text.value),
                                           int(self.slabsize_text.value),
                                           int(self.vacuum_text.value), lll_reduce=True,
                                           center_slab=self.slabcenter_checkbox.value)

            list_orientations = []
            self.full_slabs = []
            i = 0

            for slab in self.all_slabs:
                orientation = str(i) + '- ' + str(slab.miller_index[0]) + ' ' + str(slab.miller_index[1]) + ' ' \
                              + str(slab.miller_index[2])

                list_orientations.append(orientation)
                self.full_slabs.append(slab)
                i += 1

                if not self.is_structure_invertible(slab):
                    orientation = str(i) + '- ' + str(slab.miller_index[0]) + ' ' + str(slab.miller_index[1]) + ' ' \
                                  + str(slab.miller_index[2]) + ' flip'
                    list_orientations.append(orientation)

                    slab = self.flip_atoms(slab)
                    self.full_slabs.append(slab)
                    i += 1

            return self.full_slabs, len(self.full_slabs), list_orientations, self.all_slabs

    def create_slab(self, entry_id, orientation, slab):
        
        self.entry_id = entry_id

        print(type(slab))
        #import json
        #with open('temp/slab.json', 'w') as file:
         #   json.dump(slab.as_dict(), file)

        atoms = AseAtomsAdaptor.get_atoms(slab)
        atoms_slab = atoms.copy()

        structure = AseAtomsAdaptor.get_structure(atoms_slab)
        #print(type(structure))

        nx = int(self.supercellx_text.value)
        ny = int(self.supercelly_text.value)
        nz = int(self.supercellz_text.value)

        structure.make_supercell([nx, ny, nz])
        print(type(structure))

        return structure

    def adsorbate_sites(self, structure):

        print('Slab ', structure.formula)

        surface = AdsorbateSiteFinder(structure)
        self.ads_sites = surface.find_adsorption_sites()
        self.sites_dropdown.options = list(range(len(self.ads_sites['all'])))
        # self.sites_dropdown.value = self.ads_sites['all'][0]

    def add_adsorbate_ase(self, molecule=None, h=None, x=None, y=None, offset=None):
        
        if molecule is None:
            molecule = self.adatom_text.value
        
        from ase.build import add_adsorbate
        
        molecule = Atoms(molecule, positions=[(0., 0., 0.)])
        structure = self.slab_button.value[0]
        self.ads_slab = AseAtomsAdaptor.get_atoms(structure)
        cell = self.ads_slab.get_cell()
        
        if h is None:
            h = self.height_slider.value
        
        if x is None and y is None:
            x_cart = self.ads_sites['all'][self.sites_dropdown.value][0]
            y_cart = self.ads_sites['all'][self.sites_dropdown.value][1]
        else:
            x_cart = x*cell[0][0] + y*cell[1][0]
            y_cart = x*cell[0][1] + y*cell[1][1]
            
        # print(h, x, y)

        add_adsorbate(self.ads_slab, molecule, h, (x_cart, y_cart))
        # print(self.ads_slab)

        structure = AseAtomsAdaptor.get_structure(self.ads_slab)
        
        return structure

    def add_adsorbate_mol(self, molecule, slab=None, h=None, x=None, y=None, alpha=None, beta=None, gamma=None):

        from ase.build import add_adsorbate
        
        if slab is None:

            structure = self.slab_button.value[0]
            self.ads_slab = AseAtomsAdaptor.get_atoms(structure)
        else: 
            self.ads_slab = slab

        cell = self.ads_slab.get_cell()
        
        # molecule = Atoms(molecule, positions=[(0., 0., 0.)])

        # print(self.ads_sites)
        # print(self.ads_sites['all'][self.sites_dropdown.value])
        # print(x, y, h, alpha, beta, gamma)
        if alpha is None:
            molecule.rotate(self.alpha_slider.value, 'x')
        else:
            # print(type(alpha))
            alpha = float(alpha)
            # print(type(alpha))
            molecule.rotate(alpha, 'x')
            # print('rotation done')
            
        if beta is None:
            molecule.rotate(self.beta_slider.value, 'y')
        else:
            beta = float(beta)
            molecule.rotate(beta, 'y')
            
        if gamma is None:
            molecule.rotate(self.gamma_slider.value, 'z')
        else:
            gamma = float(gamma)
            molecule.rotate(gamma, 'z')
        
        if h is None:
            h = self.height_slider.value
            
        if x is None and y is None:
            x_cart = self.ads_sites['all'][self.sites_dropdown.value][0]
            y_cart = self.ads_sites['all'][self.sites_dropdown.value][1]
        else:
            x_cart = x*cell[0][0] + y*cell[1][0]
            y_cart = x*cell[0][1] + y*cell[1][1]

        add_adsorbate(self.ads_slab, molecule, h, (x_cart, y_cart))
        # print(self.ads_slab)
        
        self.ads_slab.center(axis=2, vacuum=10.0)
      
        structure = AseAtomsAdaptor.get_structure(self.ads_slab)
        
        return structure
    
    def add_adsorbate(self):

        molecule = self.adatom_text.value
        adsorbate = Molecule(molecule, [[0, 0, 0]])

        # not working, maybe use ase adsorbate instead of pymatgen
        # adsorbate_ase = Atoms(molecule, positions=[(0., 0., 0.)])
        
        # from pymatgen.io.ase import AseAtomsAdaptor
        # adsorbate = AseAtomsAdaptor.get_structure(adsorbate_ase)
        
        surface = AdsorbateSiteFinder(self.all_slabs[0])
        ads_sites = surface.find_adsorption_sites()
        
        ads_structs = surface.generate_adsorption_structures(adsorbate, repeat=[1, 1, 1])
        
        return ads_structs
