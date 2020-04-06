from ipywidgets import widgets
from IPython.display import Image, display, clear_output
from traitlets import traitlets
from ipywidgets import Layout
import numpy as np
from ase import Atoms
import json
import os

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
__date__ = "November 2019"


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class Database(object):

    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ'):
        self.rester = MPRester(MAPI_KEY)
        self.__create_showdb_button()
        self.__create_db_buttons()

        self.json_file = 'database/database.json'
        self.db = {'bulk': {}, 'slab': {}}
        if os.path.isfile('database/database.json'):
            self.__load_db("database/database.json")

    def __create_showdb_button(self):

        self.showdb_button = widgets.Button(description='Show Database',
                                            disabled=False,
                                            button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                                            tooltip='Show chosen database already calculated')

    def __create_db_buttons(self):

        self.db_buttons = widgets.ToggleButtons(options=['bulk', 'slab'],
                                                value='bulk',
                                                disabled=False)

    def __load_db(self, db_file):
        import json
        json_data = open(db_file)
        self.db = json.load(json_data)

    def load(self):
        return self.db

    def add_entry(self, mpid, energy, structure, bulk=True, slab=False):

        if bulk:
            self.db['bulk'][mpid] = energy
        elif slab:
            self.db['slab'][structure] = energy

    def write_to_json(self):

        with open(self.json_file, 'w') as outfile:
            json.dump(self.db, outfile)

    def pprint(self):
        import pprint
        with open(self.json_file, 'r') as f:
            data = f.read()
            json_data = json.loads(data)

        pprint.pprint(json_data[self.db_buttons.value])