from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Wulff construction to create the nanoparticle
"""

from six.moves import range

import itertools
from math import gcd
from functools import reduce

import numpy as np

from pymatgen.util.coord import in_coord_list

from mpinterfaces import get_struct_from_mp

from ipywidgets import widgets
from IPython.display import Image, display, clear_output
from traitlets import traitlets
from ipywidgets import Layout
import numpy as np
from ase import Atoms

from pymatgen import MPRester, Composition, Element, Molecule, Structure
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
__date__ = "April 2020"

class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))
        


class Nanoparticle_MP(Molecule):
    """
    Construct nanoparticle using wulff construction
    """

    def __init__(self, structure, rmax=15, hkl_family=[(1, 0, 0), (1, 1, 1)],
                 surface_energies=[28, 25]):

        self.structure = structure
        self.rmax = rmax
        self.hkl_family = hkl_family
        self.surface_energies = surface_energies

        spherical_neighbors = self.structure.get_sites_in_sphere(
            [0.0, 0.0, 0.0], self.rmax)
        recp_lattice = self.structure.lattice.reciprocal_lattice_crystallographic
        self.recp_lattice = recp_lattice.scale(1)
        self.set_miller_family()
        Molecule.__init__(self, [sn[0].species
                                 for sn in spherical_neighbors],
                          [sn[0].coords for sn in spherical_neighbors],
                          charge=0)

    def set_miller_family(self):
        """
        get all miller indices for the given maximum index
        get the list of indices that correspond to the given family
        of indices
        """
        recp_structure = Structure(self.recp_lattice, ["H"], [[0, 0, 0]])
        analyzer = SpacegroupAnalyzer(recp_structure, symprec=0.001)
        symm_ops = analyzer.get_symmetry_operations()

        max_index = max(max(m) for m in self.hkl_family)
        r = list(range(-max_index, max_index + 1))
        r.reverse()
        miller_indices = []
        self.all_equiv_millers = []
        self.all_surface_energies = []
        for miller in itertools.product(r, r, r):
            if any([i != 0 for i in miller]):
                d = abs(reduce(gcd, miller))
                miller_index = tuple([int(i / d) for i in miller])
                for op in symm_ops:
                    for i, u_miller in enumerate(self.hkl_family):
                        if in_coord_list(u_miller, op.operate(miller_index)):
                            self.all_equiv_millers.append(miller_index)
                            self.all_surface_energies.append(
                                self.surface_energies[i])

    def get_normals(self):
        """
        get the normal to the plane (h,k,l)
        """
        normals = []
        for hkl in self.all_equiv_millers:
            normal = self.recp_lattice.matrix[0, :] * hkl[0] + \
                self.recp_lattice.matrix[1, :] * hkl[1] + \
                self.recp_lattice.matrix[2, :] * hkl[2]
            normals.append(normal / np.linalg.norm(normal))
        return normals

    def get_centered_molecule(self):
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center
        return Molecule(self.species, new_coords,
                        charge=self._charge,
                        spin_multiplicity=self._spin_multiplicity,
                        site_properties=self.site_properties)

    def create(self):
        """
        creates the nanoparticle by chopping of the corners normal to the
        specified surfaces.
        the distance to the surface from the center of the particel =
        normalized surface energy * max radius
        """
        mol = self.get_centered_molecule()
        normalized_surface_energies = \
            np.array(self.all_surface_energies) / float(max(self.all_surface_energies))

        surface_normals = self.get_normals()
        remove_sites = []
        for i, site in enumerate(mol):
            for j, normal in enumerate(surface_normals):
                n = np.array(normal)
                n = n / np.linalg.norm(n)
                if np.dot(site.coords, n) + self.rmax * \
                        normalized_surface_energies[j] <= 0:
                    remove_sites.append(i)
                    break
        self.remove_sites(remove_sites)
        # new_sites = [site for k, site in enumerate(mol) if k not in remove_sites]
        # return Molecule.from_sites(new_sites)



    """
    Wulff construction using the ASE package
    works only for cubic systems and doesn't support multiatom basis
    from ase.cluster import wulff_construction
    from pymatgen.io.aseio import AseAtomsAdaptor
    symbol = 'Pt'
    surfaces = [ (1,0,0), (1,1,1) ]
    surface_energies = [1, 1]
    size = 200 #number of atoms
    structure = "fcc"
    latticeconstant = 5.0
    atoms = wulff_construction(symbol, surfaces, surface_energies, size, structure,
    rounding='closest', latticeconstant=latticeconstant,
    debug=False, maxiter=100)
    #convert to pymatgen structure
    pgen_structure = AseAtomsAdaptor().get_structure(atoms)
    pgen_structure.to(fmt='poscar', filename='POSCAR_pt_nano.vasp')
    """

class Nanoparticle(object):
    """
    Construct nanoparticle using wulff construction
    """
    
    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ'):

        self.rester = MPRester(MAPI_KEY)

        self.__create_rmax_text()
        self.__create_surface_families_text()
        self.__create_surface_energies_text()

        self.__create_nanoparticle_button()


        #self.structure = structure
        self.rmax = None
        self.list_families = None
        self.list_families =None
        

    def __create_rmax_text(self):
        """Build the Text widget for Miller index  input"""
        self.rmax_text = widgets.Text(    
            placeholder='rmax',
            value='25',
            description='rmax',
            disabled=False,
            layout=Layout(width='150px'),
            style={'description_width': 'initial'}
        )

    def __create_surface_families_text(self):
        """Build the Text widget for Miller index  input"""
        self.surface_families_text = widgets.Text(    
            placeholder='(0, 0, 1) ; (0, 1, 1)',
            value='(0, 0, 1) ; (0, 1, 1)',
            description='Surface families',
            disabled=False,
            layout=Layout(width='500px'),
            style={'description_width': 'initial'}
        )

    def read_surface_families(self):
        list_families = self.surface_families_text.value.split(';')
        self.list_families = []

        for family in list_families:
            self.family = []
            family = family.replace(" ", "")
            family = family.replace("(", "")
            family = family.replace(")", "")
            list_index = family.split(',')

            for index in list_index:
                self.family.append(int(index))

            self.list_families.append(tuple(self.family))

    def read_rmax(self):
        self.rmax = int(self.rmax_text.value)

    def read_surface_energies(self):
        list_energies = self.surface_energies_text.value.split(',')
        self.list_energies = []

        for energy in list_energies:
            
            energy = energy.replace(" ", "")
            self.list_energies.append(int(energy))

    def __create_surface_energies_text(self):
        """Build the Text widget for Miller index  input"""
        self.surface_energies_text = widgets.Text(    
            placeholder='18, 20',
            value='18, 20',
            description='Surface energies',
            disabled=False,
            layout=Layout(width='500px'),
            style={'description_width': 'initial'}
        )


    def __create_nanoparticle_button(self):
        """Build the button widget to create nanoparticle."""

        # create widget
        self.nanoparticle_button = LoadedButton(description='Create Nanoparticle', value=None, disabled=False)


    def create_nanoparticle(self, mp_id):

        self.read_surface_families()
        self.read_surface_energies()
        self.read_rmax()

        structure = self.rester.get_structure_by_material_id(mp_id)

        sa = SpacegroupAnalyzer(structure)
        structure_conventional = sa.get_conventional_standard_structure()

        self.nanop = Nanoparticle_MP(structure_conventional, rmax=self.rmax,
                                        hkl_family=self.list_families,
                                        surface_energies=self.list_energies)

        self.nanop.create()

        self.nanop.to(fmt='xyz', filename='output_cif/nanoparticle.xyz')

        from ase.io import read, write

        np_atoms = read('output_cif/nanoparticle.xyz')

        return np_atoms


    