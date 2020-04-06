from ipywidgets import widgets
from IPython.display import Image, display, clear_output
from traitlets import traitlets

from pymatgen import MPRester, Composition, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility


compatibility = MaterialsProjectCompatibility()

import json
import os
import sys

if sys.version_info[0] == 2:
    from urllib import quote_plus
else:
    from urllib.parse import quote_plus

import requests

urlpattern = {
    "results": "https://materialsproject.org/molecules/results?query={spec}",
    "mol_json": "https://materialsproject.org/molecules/{mol_id}/json",
    "mol_svg": "https://materialsproject.org/molecules/{mol_id}/svg",
    "mol_xyz": "https://materialsproject.org/molecules/{mol_id}/xyz",
}


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


class Search(object):
    
    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ'):
        
        """Initialize elements of the GUI.
           Arguments
           ---------

        """

        self.rester = MPRester(MAPI_KEY)
        
        # Create __structure Input Text
        self.__create_structure()
        self.__create_database_dropdown()
        self.__create_composition_checkbox()
        self.__create_molecule_checkbox()
        self.__create_search_button()
        self.__create_search_select()
        self.__create_unitcell_radio()
        self.__create_viewer_radio()
        self.__create_viewermol_radio()
        self.__create_phasediagram_button()
        self.__create_mapikey_text()
        
    def __create_structure(self):
        """Build the Text widget for composition input"""
        self.structure = widgets.Text(    
            placeholder='Type something',
            value='Al2Cu',
            description='Structure',
            disabled=False,
            style={'description_width': 'initial'}
        )

    def __create_mapikey_text(self):
        self.mapi_key = widgets.Text(
            placeholder='Materials Project API Key',
            value='4oBTKz0pkFSg9EUQ',
            description='MAPI Key',
            disabled=False,
            style={'description_width': 'initial'}
        )

    def __create_database_dropdown(self):
        """Build the Dropdowm widget for database choice"""
        self.database_dropdown = widgets.Dropdown(options=['Materials Project'],
                                                  value='Materials Project',
                                                  description = 'Database',
                                                  style={'description_width': 'initial'})

    def __create_composition_checkbox(self):
        """Build the Checkbox widget for Exact composition search."""
        self.composition_checkbox = widgets.Checkbox(value = False, 
                                                     description='Chemical Search',
                                                     style={'description_width': 'initial'})

    def __create_molecule_checkbox(self):
        """Build the Checkbox widget for molecule search."""
        self.molecule_checkbox = widgets.Checkbox(value = False, 
                                                  description = 'Molecule',
                                                  style={'description_width': 'initial'})

    def __create_unitcell_radio(self):
        """Build the Radio widget for Unit Cell choice."""
        self.unitcell_radio = widgets.RadioButtons(options=['primitive', 'conventional'], 
                                                   value='primitive', 
                                                   description='Unit Cell:',
                                                   disabled=False)
        
    def __create_viewer_radio(self):
        """Build the Radio widget for display outputs."""
        self.viewer_radio = widgets.RadioButtons(options=['3D structure', 'Bands structure', 'DOS'],
                                                 value=None,
                                                 description='Output:',
                                                 disabled=False)    

    def __create_viewermol_radio(self):
        """Build the Radio widget for display outputs."""
        self.viewermol_radio = widgets.RadioButtons(options=['3D structure'],
                                                    value='3D structure',
                                                    description='Output:',
                                                    disabled=True)
        
    def __create_search_button(self):
        """Build the button widget to search materials."""
        
        # create widget
        self.search_button = LoadedButton(description='Search', value=False,
                                          tooltip='Search the material on the selected database', button_style='info')
        
        # update function
        self.search_button.on_click(self.search_clicked)
        
    def __create_search_select(self):
        """Build the Select widget for search results ."""
        
        # create widget
        self.search_select = widgets.Select(options=['None'],
                                            value='None',
                                            rows=5,
                                            description='')

    def __create_phasediagram_button(self):
        """Build the button widget for Phase Diagram Generation."""
        
        # create widget
        self.phasediagram_button = LoadedButton(description='Phase Diagram', value=None, disabled=True,
                                                button_style='danger')
        
        # update function
        # self.phasediagram_button.on_click(self.pd_clicked)

    def search_molecules(self, composition):
        # self.composition_checkbox.value = False
       
        spec = {"formula": composition}
        # Stringify `spec`, ensure the string uses double quotes, and percent-encode it...
        str_spec = quote_plus(str(spec).replace("'", '"'))
        # ...because the spec is the value of a "query" key in the final URL.
        url = urlpattern["results"].format(spec=str_spec)
            
        return (requests.get(url, headers={'X-API-KEY': MAPI_KEY})).json()

    def get_molecule(self, mol_id, fmt='json'):
        url = urlpattern["mol_" + fmt].format(mol_id=mol_id)
        response = requests.get(url, headers={'X-API-KEY': MAPI_KEY})
        if fmt == 'json':
            return response.json()
        else:
            return response.content 
    
    def search_clicked(self, b):
        b.value = None
        # print('Hello')
        composition = self.structure.value
        
        if self.molecule_checkbox.value:
            
            self.viewer_radio.options = ['3D structure']
            self.viewer_radio.value = '3D structure'
            self.viewer_radio.disabled = True
            
            # Query for a molecule
            # Call Search_Molecule fonction 
            print('Not done yet')
            list_molecules = []
            list_energy = []
            dropdown = []
            results = self.search_molecules(composition)
            id_stable = results[0]['task_id']
            dropdown_stable = results[0]['task_id'] + ' {1} ({0:.4f})'.format(results[0]['IE'], results[0]['formula'])
            for molecule in results:
                list_molecules.append(molecule['task_id'])
                list_energy.append(molecule['IE'])
                dropdown.append(molecule['task_id'] + ' {1} ({0:.4f})'.format(molecule['IE'], molecule['formula']))

            if id_stable is None:
                # print('No stable structure')
                print('No molecule')
                b.value = None
                self.phasediagram_button.disabled = True
                self.phasediagram_button.button_style = 'danger'
            else:
                b.value = [id_stable, results, list_molecules, list_energy]
                
                self.phasediagram_button.disabled = True
                self.phasediagram_button.button_style = 'danger'
               
                self.search_select.options = dropdown
                self.search_select.value = dropdown_stable
        else:
            
            # Query for a crystal structure
            self.viewer_radio.options = ['3D structure', 'Bands structure', 'DOS']
            self.viewer_radio.value = None
            self.viewer_radio.disabled = False

            if self.composition_checkbox.value is False:
                # Exact Composition
                
                list_structures = []
                list_ehull = []
                dropdown = []
                
                # Import all entries of MP with composition
                try:
                    mp_entries = self.rester.get_entries(composition, sort_by_e_above_hull=True,
                                                    property_data=['pretty_formula'])
                except:
                    print("No acces to Materials Project Database")
                
                # Search the stable structure 
                
                min_ehull = mp_entries[0].data['e_above_hull']
                id_stable = mp_entries[0].entry_id
                dropdown_stable = mp_entries[0].entry_id + ' {1} ({0:.4f})'.format(mp_entries[0].data['e_above_hull'],
                                                                                   mp_entries[0].data['pretty_formula'])
                
                for entry in mp_entries:
                    list_structures.append(entry.entry_id)
                    list_ehull.append(entry.data['e_above_hull'])
                    dropdown.append(entry.entry_id + ' {1} ({0:.4f})'.format(entry.data['e_above_hull'],
                                                                             entry.data['pretty_formula']))

                    if entry.data['e_above_hull'] < min_ehull:
                        id_stable = entry.entry_id
                        dropdown_stable = entry.entry_id + ' {1} ({0:.4f})'.format(entry.data['e_above_hull'],
                                                                                   entry.data['pretty_formula'])
                        
                if id_stable is None:
                    print('No stable structure')
                    b.value = None
                    
                else:
                    b.value = [id_stable, mp_entries, list_structures, list_ehull]
                    self.search_select.options = dropdown
                    self.search_select.value = dropdown_stable
                    self.phasediagram_button.disabled = True
                    self.phasediagram_button.button_style = 'danger'
                
            else:
                # Chemical Elements
                compo = Composition(composition)
                list_elements = []
                for i in compo.elements:
                    list_elements.append(i.value)
                # Get all entries in the chemical system
                mp_entries = self.rester.get_entries_in_chemsys(list_elements, property_data=['pretty_formula'])
                entries = compatibility.process_entries(mp_entries)

                list_structures = []
                list_ehull = []
                dropdown = []
                list_stab = self.rester.get_stability(entries)
                min_ehull = list_stab[0]['e_above_hull']
                id_stable = mp_entries[0].entry_id
                dropdown_stable = mp_entries[0].entry_id + ' {1} ({0:.4f})'.format(list_stab[0]['e_above_hull'],
                                                                                   mp_entries[0].data['pretty_formula'])

                for i, entry in enumerate(mp_entries):
                    list_structures.append(entry.entry_id)
                    list_ehull.append(list_stab[i]['e_above_hull'])
                    dropdown.append(entry.entry_id + ' {1} ({0:.4f})'.format(list_stab[i]['e_above_hull'],
                                                                             entry.data['pretty_formula']))
                    # dropdown.append(structure['entry_id'])

                if id_stable is None:
                    print('No stable structure')
                    b.value = None
                    
                else:
                    b.value = [id_stable, mp_entries, list_structures, list_ehull]
                    self.search_select.options = dropdown
                    self.search_select.value = dropdown_stable
                    self.phasediagram_button.disabled = False
                    self.phasediagram_button.style.button_color = 'lightgreen'
