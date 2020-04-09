import sys
import json

try:
    import GPyOpt
    from GPyOpt.methods import BayesianOptimization
except:
    sys.exit("GPyOpt is necessary. Run pip install gpyopt")
    
try:
    import ase
    from ase import Atoms
    #from ase.build import fcc111, add_adsorbate, bulk
    from ase.io import read, write
except:
    sys.exit("ase is necessary. Run pip install --upgrade --user ase")
    
try:
    from paramiko import SSHClient
    import paramiko
except:
    sys.exit("paramiko is necessary. Run pip install paramiko")

    
try:
    from pymatgen import MPRester, Composition, Element
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.electronic_structure.plotter import BSPlotter
    from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, CompoundPhaseDiagram
    from pymatgen.analysis.phase_diagram import PDPlotter
    from pymatgen.io.vasp import Vasprun
    from pymatgen.entries.computed_entries import ComputedEntry
    from pymatgen.entries.compatibility import MaterialsProjectCompatibility
    from pymatgen.util.plotting import pretty_plot
    from pymatgen.analysis.adsorption import *
except:
    sys.exit("pymatgen is necessary. Run conda install --channel "
             "matsci pymatgen and conda upgrade pymatgen Or pip install --upgrade pymatgen")
    
import numpy as np
import ipywidgets
from ipywidgets import widgets
from IPython.display import display, clear_output
from ipywidgets import Layout
from traitlets import traitlets
from matplotlib import pyplot as plt
from ase.visualize import view
from ipyfilechooser import FileChooser

viewer = 'ngl'

import importlib
from src import optimizer_module, remote_module, viewer_mod, search_module, builder_module, vasp_module, data_module, qe_module, np_module, aimd_module, optmd_module

# Change in mymodule
importlib.reload(search_module)
importlib.reload(builder_module)
importlib.reload(optimizer_module)
importlib.reload(remote_module)
importlib.reload(vasp_module)
importlib.reload(data_module)
importlib.reload(viewer_mod)
importlib.reload(qe_module)
importlib.reload(np_module)
importlib.reload(aimd_module)
importlib.reload(optmd_module)


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


class GUI(object):
      
    def __init__(self, MAPI_KEY='4oBTKz0pkFSg9EUQ', VASP_PP_PATH="/home/fgimbert/Projects/potvasp", ESPRESSO_PSEUDO="/home/fgimbert/Projects/potvasp"):
        """Initialize elements of the GUI.
           Arguments
           ---------

        """
        self.MAPI_KEY = MAPI_KEY
        self.VASP_PP_PATH = VASP_PP_PATH
        self.ESPRESSO_PSEUDO = ESPRESSO_PSEUDO
        self.rester = MPRester(self.MAPI_KEY)

        self.emin = None
        self.slabmin = None
        self.domain = []

        self.mpid = None # The mpid for the surface material - not used now
        self.current_structure = None # The structure (bulk or slab) displayed
        self.all_orientations = None # All the orientations for a slab depending of the input Miller index
        self.all_slabs = None # All the slabs
        self.current_slab = None  # The slab generated for a specific orientation chosen in dropdown
        self.orientation_ckecked = False

        self.remote_workdir = None

        self.search_menu = search_module.Search(MAPI_KEY=self.MAPI_KEY)
        self.builder = builder_module.Builder(MAPI_KEY=self.MAPI_KEY)
        self.nanoparticle = np_module.Nanoparticle(MAPI_KEY=self.MAPI_KEY)
        self.aimd = aimd_module.AIMD(MAPI_KEY=self.MAPI_KEY, VASP_PP_PATH=self.VASP_PP_PATH, ESPRESSO_PSEUDO=self.ESPRESSO_PSEUDO)
        self.optmd = optmd_module.OptMD(MAPI_KEY=self.MAPI_KEY, VASP_PP_PATH=self.VASP_PP_PATH, ESPRESSO_PSEUDO=self.ESPRESSO_PSEUDO)

        self.optimizer = optimizer_module.Optimizer(MAPI_KEY=self.MAPI_KEY)
        self.remote = remote_module.Remote()
        self.vaspconfig = vasp_module.configVASP(MAPI_KEY=self.MAPI_KEY, VASP_PP_PATH=self.VASP_PP_PATH)
        self.qeconfig = qe_module.configQE(MAPI_KEY=self.MAPI_KEY, ESPRESSO_PSEUDO=self.ESPRESSO_PSEUDO)

        self.database = data_module.Database(MAPI_KEY=self.MAPI_KEY)

        self.bulk_db = self.database.load()

        self.__create_display_button()
        self.__create_savecif_button()
        self.__create_selectcif_button()
        self.__create_uploadcif_button()
        self.__create_saveposcar_button()
        self.__create_unitcellx_slider()
        self.__create_unitcelly_slider()
        self.__create_unitcellz_slider()
        self.__create_calcbuttons()

        self.search_menu.phasediagram_button.on_click(self.pd_clicked)

        self.builder.orientation_button.on_click(self.checkorien_clicked)
        self.builder.slab_button.on_click(self.slab_clicked)
        self.builder.sernergy_button.on_click(self.surfacenergy_clicked)
        self.builder.adsorbate_button.on_click(self.add_adsorbate_clicked)
        self.builder.molad_button.on_click(self.checkmol_clicked)

        self.nanoparticle.nanoparticle_button.on_click(self.nanoparticle_clicked)

        self.aimd.rundft_button.on_click(self.rundft_clicked)
        self.optmd.runoptmd_button.on_click(self.runoptmd_clicked)


        self.optimizer.runopt_button.on_click(self.runopt_clicked)
        # self.optimizer.testopt_button.on_click(self.testopt_clicked)

        self.optimizer.runcalc_button.on_click(self.runcalc_clicked)
        self.optimizer.inputvasp_button.on_click(self.createinput_clicked)

        self.database.showdb_button.on_click(self.showdb_clicked)

        self.remote.jobsrunning_button.on_click(self.jobsrunning_clicked)
        self.remote.boupdate_button.on_click(self.boupdate_clicked)
        self.remote.outputjob_button.on_click(self.booutput_clicked)
        self.remote.killjob_button.on_click(self.killjob_clicked)

        self.input_panel = ipywidgets.Output(layout=Layout(width='100%', border='solid 1px'))
        self.selection_panel = ipywidgets.Output(layout=Layout(width='50%', border='solid 1px', height='230px'))
        self.description_panel = ipywidgets.Output(layout=Layout(width='50%', border='solid 1px',
                                                                 height='230px', overflow_y='auto'))


        self.import_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                            height='250px', overflow_y='auto'))

        self.building_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                              height='320px', overflow_y='auto'))

        self.nanopart_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                              height='320px', overflow_y='auto'))

        self.aimd_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                              height='320px', overflow_y='auto'))  
        self.optmd_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                              height='320px', overflow_y='auto')) 

        self.slab_panel = ipywidgets.Output(layout=Layout(width='60%', height='270px', overflow_y='auto'))
        self.sites_panel = ipywidgets.Output(layout=Layout(width='40%', border='solid 1px',
                                                           height='300px', overflow_y='auto'))

        self.adsorbate_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                               height='320px', overflow_y='auto'))
        self.optimization_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                                  height='320px', overflow_y='auto'))

        self.config_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                            height='750px', overflow_y='auto'))
        self.jobs_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                          height='320px', overflow_y='auto'))
        self.jobslist_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 0px', height='260px',
                                                              overflow_y='auto'))

        self.data_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                          height='320px', overflow_y='auto'))

        self.showdb_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px', height='180px',
                                                            overflow_y='auto'))

        self.remote_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px', height='200px',
                                                            overflow_y='auto'))

        self.vaspconfig_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px', height='200px',
                                                                overflow_y='auto'))

        self.qeconfig_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px', height='200px',
                                                                overflow_y='auto'))

        self.jobconfig_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px', height='300px',
                                                                overflow_y='auto'))

        self.inputopt_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 0.5px',
                                                              height='190px', overflow_y='auto'))

        self.outputopt_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 0.5px',
                                                               height='410px', overflow_y='auto'))

        self.outputopt2_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 0.5px',
                                                                height='410px', overflow_y='auto'))

        self.accordion_outopt = widgets.Accordion(children=[self.outputopt_panel,
                                                            self.outputopt2_panel], selected_index=None)
        self.accordion_outopt.set_title(0, 'BO Run Steps')
        self.accordion_outopt.set_title(1, 'BO Output')

        self.output_panel = ipywidgets.Output(layout=Layout(width='80%', border='solid 2px',
                                                            height='500px', overflow_y='auto'))
        self.action_panel = ipywidgets.Output(layout=Layout(width='20%', border='solid 2px',
                                                            height='500px', overflow_y='auto'))
        
        self.outputs_panels = widgets.HBox([self.action_panel, self.output_panel])

        self.ssh_panel = ipywidgets.Output(layout=Layout(width='99%', border='solid 1px', height='500px',
                                                         overflow_y='auto'))
        
        self.error_panel = ipywidgets.Output(layout=Layout(width='99%', border='solid 1px', height='100px',
                                                           overflow_y='auto'))
        self.molecules_panel = ipywidgets.Output(layout=Layout(width='60%', height='270px', overflow_y='auto'))
        self.adsorb_panel = ipywidgets.Output(layout=Layout(width='40%', border='solid 1px', height='300px',
                                                            overflow_y='auto'))
        
        self.building_tab = widgets.VBox([self.building_panel])
        self.adsorbate_tab = widgets.VBox([self.adsorbate_panel])
        self.optimizer_tab = widgets.VBox([self.optimization_panel])
        self.config_tab = widgets.VBox([self.config_panel])
        self.jobs_tab = widgets.VBox([self.jobs_panel])
        self.data_tab = widgets.VBox([self.data_panel])

        self.structure_panel = widgets.HBox([self.selection_panel, self.description_panel])
        self.input_tab = widgets.VBox([self.input_panel,self.structure_panel])
        self.import_tab = widgets.VBox([self.import_panel])

        

        # Tabs for Sla optimization

        tab_contents = [self.building_tab, self.adsorbate_tab, self.optimizer_tab]
        tab_name = ['1. Create Slab', '2. Add Adsorbate', '3. Optimize Structure']

        children = [content for content in tab_contents]

        self.tab = widgets.Tab()

        self.tab.children = children
        for i in range(len(children)):
            self.tab.set_title(i, tab_name[i])


        # Tabs for Materials Informatics Tools

        self.slabopt_tab = widgets.VBox([self.tab])
        self.aimd_tab = widgets.VBox([self.aimd_panel])
        self.optmd_tab = widgets.VBox([self.optmd_panel])
        self.nanopart_tab = widgets.VBox([self.nanopart_panel])

        miprocess_tabs = [self.slabopt_tab, self.nanopart_tab, self.aimd_tab, self.optmd_tab, self.jobs_tab, self.data_tab]
        miprocess_name = ['Slab Optimization', 'Nanoparticles', 'AIMD', 'On-the-fly MD', 'Current Jobs', 'Database']
        children = [content for content in miprocess_tabs]

        self.miprocess_tab = widgets.Tab()
        self.miprocess_tab.children = children
        for i in range(len(children)):
            self.miprocess_tab.set_title(i, miprocess_name[i])

        menu_tab = [self.input_tab, self.import_tab, self.config_tab]
        menu_name = ['Search Materials', "Import Structure", 'Configuration']
        self.menu = widgets.Tab()
        menu_children = [content for content in menu_tab]
        self.menu.children = menu_children
        for i in range(len(menu_children)):
            self.menu.set_title(i, menu_name[i])

        text_search = '1. Creation Materials Tools:'
        searchWidget = widgets.HTML(value=f"<b><font color='red'>{text_search}</b>")

        self.display_input_panel()
        self.display_import_panel()
        self.display_action_panel()
            
        with self.selection_panel:
            text_search = 'Search Results:'
            text_pd = 'If chemical search:'
            resultsWidget = widgets.HTML(value=f"<b><font color='black'>{text_search}</b>")
            pdWidget = widgets.HTML(value=f"<b><font color='black'>{text_pd}</b>")
         
            results_box = widgets.VBox([resultsWidget,
                                        widgets.HBox([self.search_menu.search_select, self.search_menu.unitcell_radio]),
                                        widgets.HBox([self.search_menu.viewer_radio,
                                                      widgets.VBox([self.__display_button, pdWidget,
                                                                    self.search_menu.phasediagram_button])])])

            display(results_box)
        
        text_building = '2. Materials Informatics Tools:'
        buildingWidget = widgets.HTML(value=f"<b><font color='red'>{text_building}</b>")
        
        self.display_building_panel()
        self.display_adsorbate_panel()

        self.display_nanopart_panel()
        self.display_aimd_panel()
        self.display_optmd_panel()
        self.display_optimizer_panel()

        self.display_config_panel()
        self.display_data_panel()
        self.display_jobs_panel()

        # self.display_inputopt_panel()
        
        text_optimization = '3. Optimization Structure:'
        optimizationWidget = widgets.HTML(value=f"<b><font color='red'>{text_optimization}</b>")
        with self.optimization_panel:
            optimization_box = widgets.VBox([optimizationWidget])
            # display(optimization_box)
        
        text_output = '3D Display:'
        outputWidget = widgets.HTML(value=f"<b><font color='red'>{text_output}</b>")
        
        with self.output_panel:
            results_box = widgets.VBox([outputWidget])
            # display(results_box)
        
        text_error = 'Error Output:'
        errorWidget = widgets.HTML(value=f"<b><font color='red'>{text_error}</b>")

        text_ssh = 'Output'
        sshWidget = widgets.HTML(value=f"<b><font color='red'>{text_ssh}</b>")

        with self.error_panel:
            # print(self.remote.host)
            self.remote.read_remotehosts()
            error_box = widgets.VBox([errorWidget])
            # display(error_box)
            
        self.search_menu.search_select.observe(self.update_description_panel, 'value')
        self.search_menu.molecule_checkbox.observe(self.update_slab_button, 'value')

        # self.builder.orientation_select.observe(self.update_slab_button, 'value')
        self.search_menu.search_select.observe(self.update_orientation_button, 'value')

        self.search_menu.unitcell_radio.observe(self.update_description_panel, 'value')

        self.remote.cluster_button.on_click(self.update_ssh_panel)

        # self.search_menu.search_select.observe(self.update_output_panel, 'value')
        # self.search_menu.unitcell_radio.observe(self.update_output_panel, 'value')
        # self.search_menu.viewer_radio.observe(self.update_output_panel, 'value')
        self.search_menu.search_button.observe(self.update_display_button, 'value')  
        self.search_menu.viewer_radio.observe(self.update_display_button, 'value') 

        self.__uploadcif_button.observe(self.update_display_button, 'value')

        self.__display_button.observe(self.update_output_panel, 'value')
        # self.__display_button.observe(self.update_sites_panel, 'value')
        self.__display_button.observe(self.update_unitcell_sliders, 'value')

        self.builder.adatom_text.observe(self.update_adsorbate_button, 'value')
        # self.search_menu.phasediagram_button.observe(self.update_output_panel, 'value')



        display(widgets.VBox([searchWidget,
                              self.menu,
                              buildingWidget, 
                              self.miprocess_tab,
                              # optimizationWidget,
                              # self.optimization_panel,
                              outputWidget, 
                              self.outputs_panels,
                              sshWidget, self.ssh_panel,
                              errorWidget, self.error_panel]))
    
    def __create_savecif_button(self):
        """Build the Save CIF button widget."""
        
        # create widget
        self.__savecif_button = LoadedButton(description='Save CIF', value=None, disabled=True, button_style='warning')
        
        # update function
        self.__savecif_button.on_click(self.savecif_clicked)

    def __create_selectcif_button(self):
        """Build the Upload CIF button widget."""
        self.__selectcif_button = FileChooser()
        # update function
        #self.uploadcif_button.on_click(self.uploadcif_clicked)
    
    def __create_uploadcif_button(self):
        """Build the Upload CIF button widget."""


        self.__uploadcif_button = LoadedButton(description='Display Structure', value=None, disabled=False, button_style='warning')
        
        # update function
        self.__uploadcif_button.on_click(self.uploadcif_clicked)


    def __create_saveposcar_button(self):

        """Build the Save CIF button widget."""

        # create widget
        self.__saveposcar_button = LoadedButton(description='Save POSCAR', value=None, disabled=True,
                                                button_style='warning')

        # update function
        self.__saveposcar_button.on_click(self.saveposcar_clicked)


    def __create_unitcellx_slider(self):
        """Build the unit cell x slider widget."""
        self.unitcellx_slider = widgets.IntSlider(value=1,
                                                  min=1,
                                                  max=15,
                                                  step=1,
                                                  description='x:',
                                                  layout=Layout(width='170px'),
                                                  style={'description_width': 'initial'},
                                                  disabled=True
                                                  )

    def __create_unitcelly_slider(self):
        """Build the unit cell x slider widget."""
        self.unitcelly_slider = widgets.IntSlider(value=1,
                                                  min=1,
                                                  max=15,
                                                  step=1,
                                                  description='y:',
                                                  layout=Layout(width='170px'),
                                                  style={'description_width': 'initial'},
                                                  disabled=True
                                                  )

    def __create_unitcellz_slider(self):
        """Build the unit cell x slider widget."""
        self.unitcellz_slider = widgets.IntSlider(value=1,
                                                  min=1,
                                                  max=15,
                                                  step=1,
                                                  description='z:',
                                                  layout=Layout(width='170px'),
                                                  style={'description_width': 'initial'},
                                                  disabled=True
                                                  )

    def __create_calcbuttons(self):
        self.calc_buttons = widgets.ToggleButtons(options=['EMT', 'VASP', 'QE (in work)'],
                                                  value='EMT',
                                                  tooltips=['Python calculator',
                                                            'VASP calculator. Needs cluster configuration',
                                                            'QE calculator. In work'])

    def savecif_clicked(self, b):
        
        from ase.io import write
        name = 'slab'
        write('output_cif/{0}.cif'.format(name), self.__savecif_button.value)

    def uploadcif_clicked(self, b):

        self.ciffile = self.__selectcif_button.selected

        self.__display_button.value = ['Input file',  self.ciffile]
        
        with self.error_panel:
            clear_output()
            print('Button Upload Clicked !')
            print(self.__display_button.value)
        

    def saveposcar_clicked(self, b):

        from ase.io import write
        name = 'POSCAR'
        write('output_cif/{0}'.format(name), self.__savecif_button.value, format='vasp')

    def update_unitcell_sliders(self, *args):

        if self.__display_button.value[1] == '3D structure':
            self.unitcellx_slider.disabled = False
            self.unitcellz_slider.disabled = False
            self.unitcelly_slider.disabled = False
        else:
            self.unitcellx_slider.disabled = True
            self.unitcellz_slider.disabled = True
            self.unitcelly_slider.disabled = True

    def update_ssh_panel(self, *args):

        with self.ssh_panel:
            clear_output()
            # print('Passphrase: ', self.remote.passphrase)
            # print('Password: ', self.remote.password)
            print('Cluster: ', self.remote.hostname)
            print('User: ', self.remote.user)
            # print('Port: ', self.remote.port)
            print('Host: ', self.remote.host)

            ls_success = self.remote.command_ls()
            if ls_success == 0:
                self.remote.cluster_button.icon = 'check'
                self.remote.cluster_button.description='Connected'

            # self.remote.put_file()
            # self.remote.command_ls()

    def display_action_panel(self):

        
        text = 'If bulk:'
        # input_label = widgets.Label(value = r'\(\color{red} {' + text  + '}\)')
        bulkWidget = widgets.HTML(value=f"<b><font color='black'>{text}</b>")
        """Display the input panel for Search Module"""
        with self.action_panel:
            input_search = widgets.VBox([bulkWidget,
                                        self.unitcellx_slider,
                                        self.unitcelly_slider,
                                        self.unitcellz_slider,
                                        self.__savecif_button,
                                        self.__saveposcar_button
                                        ])

            # input_box = widgets.VBox([htmlWidget, input_search])
            display(input_search)

    def display_input_panel(self):
        
        """Display the input panel for Search Module"""
        with self.input_panel:
            text = '1. Search Material:'
            # input_label = widgets.Label(value = r'\(\color{red} {' + text  + '}\)')
            htmlWidget = widgets.HTML(value=f"<b><font color='red'>{text}</b>")
            input_search = widgets.VBox([widgets.HBox([self.search_menu.structure,
                                         self.search_menu.database_dropdown,
                                         self.search_menu.composition_checkbox,
                                         self.search_menu.molecule_checkbox,
                                         self.search_menu.search_button])])
            
            input_box = widgets.VBox([input_search])
            
            display(input_box)
            
    def display_building_panel(self):
        
        """Display the building panel for Builder Module"""
        
        with self.building_panel:
            
            with self.slab_panel:

                text_slab = 'Slab parameters:'
                text_cell = 'Supercell size:'
                text_miller = 'Slab orientation'
                # input_label = widgets.Label(value = r'\(\color{red} {' + text  + '}\)')
                slabWidget = widgets.HTML(value=f"<b><font color='black'>{text_slab}</b>")
                cellWidget = widgets.HTML(value=f"<b><font color='black'>{text_cell}</b>")
                millerWidget = widgets.HTML(value=f"<b><font color='black'>{text_miller}</b>")
            
                input_supercell = widgets.HBox([cellWidget,widgets.HBox([self.builder.supercellx_text,
                                                                         self.builder.supercelly_text,
                                                                         self.builder.supercellz_text])])
            
                input_search = widgets.VBox([widgets.HBox([self.builder.vacuum_text,self.builder.slabcenter_checkbox]),
                                             self.builder.slabsize_text])

                input_miller = widgets.HBox([self.builder.miller_text, self.builder.orientation_button,
                                             self.builder.orientation_select])

                output_box = widgets.HBox([  # self.builder.slabs_dropdown,
                                           self.builder.slab_button,
                                           self.builder.sernergy_button, self.builder.all_energies_button])

                input_box = widgets.VBox([slabWidget,
                                          input_search,
                                          input_supercell,
                                          millerWidget,
                                          input_miller,
                                          output_box])

                display(widgets.HBox([input_box]))

            with self.sites_panel:
                print("Relevant adsorption sites: ")

            display(widgets.HBox([self.slab_panel, self.sites_panel]))

    def display_adsorbate_panel(self):
        """Display the adsorbate panel for Adsorbate Module"""
        
        with self.adsorbate_panel:

            with self.molecules_panel:

                text_adatom = 'Adsorbate choice:'
                # input_label = widgets.Label(value = r'\(\color{red} {' + text  + '}\)')
                adatomText = widgets.HTML(value=f"<b><font color='black'>{text_adatom}</b>")
                
                atomWidget = widgets.HBox([self.builder.adatom_radio, 
                                           self.builder.adatom_text, self.builder.molad_button])

                display(widgets.VBox([adatomText, 
                                      atomWidget,
                                      self.builder.sites_dropdown,
                                      widgets.HBox([self.builder.height_slider,
                                                    widgets.VBox([ self.builder.alpha_slider,
                                                                   self.builder.beta_slider,
                                                                   self.builder.gamma_slider])]),
                                      self.builder.adsorbate_button]))

            with self.adsorb_panel:
                
                text_mol = 'Only if adsorbate is molecule:'
                # input_label = widgets.Label(value = r'\(\color{red} {' + text  + '}\)')
                self.molText = widgets.HTML(value=f"<b><font color='black'>{text_mol}</b>")
                
                display(widgets.VBox([self.molText,
                                      self.builder.molad_dropdown]))

            display(widgets.HBox([self.molecules_panel, self.adsorb_panel]))

    def display_nanopart_panel(self):

        with self.nanopart_panel:
            text_aimd = 'Module for nanoparticles creations (in work)'
            nanopartWidget = widgets.HTML(value=f"<b><font color='red'>{text_aimd}</b>")

            input_nanoparticle = widgets.VBox([self.nanoparticle.rmax_text, 
                                                self.nanoparticle.surface_families_text, 
                                                self.nanoparticle.surface_energies_text,
                                                self.nanoparticle.nanoparticle_button])

            input_box = widgets.VBox([nanopartWidget,
                                        input_nanoparticle])

            display(input_box)
            
            #slabWidgets = widgets.VBox([self.tab])
            #display(aimdWidget)

    def display_aimd_panel(self):

        with self.aimd_panel:
            text_aimd = 'Module for single DFT calculations (in work)'
            aimdWidget = widgets.HTML(value=f"<b><font color='red'>{text_aimd}</b>")

            input_aimd = widgets.VBox([self.aimd.vasppp_text, 
                                                self.aimd.vasprun_text,
                                                self.aimd.workdir_text,
                                                self.aimd.ncpus_text,
                                                self.aimd.rundft_button
                                                ])



            input_box = widgets.VBox([aimdWidget,
                                        input_aimd])
            #slabWidgets = widgets.VBox([self.tab])
            display(input_box)

    def display_optmd_panel(self):

        with self.optmd_panel:
            text_optmd = 'Module for on-the-fly optimized AIMD calculations (in work)'
            optmdWidget = widgets.HTML(value=f"<b><font color='red'>{text_optmd}</b>")

            input_optmd = widgets.VBox([self.optmd.vasppp_text, 
                                                self.optmd.vasprun_text,
                                                self.optmd.workdir_text,
                                                self.optmd.ncpus_text,
                                                self.optmd.runoptmd_button
                                                ])



            input_box = widgets.VBox([optmdWidget,
                                        input_optmd])
            #slabWidgets = widgets.VBox([self.tab])
            display(input_box)

    def display_optimizer_panel(self):
        """Display the optimizer panel for Optimizer Module"""
        
        with self.optimization_panel:
            
            with self.inputopt_panel:

                calcWidgets = widgets.VBox([# self.optimizer.inputvasp_button,
                                            self.optimizer.runcalc_button,
                                            self.optimizer.runopt_button,
                                            #self.optimizer.testopt_button
                                            ])
                calc2Widgets = widgets.VBox([self.optimizer.initbo_text,
                                             self.optimizer.initbo_slider,
                                             self.optimizer.maxbo_text,
                                             self.optimizer.maxbo_slider])

                display(widgets.HBox([self.optimizer.parameters_select, calc2Widgets, calcWidgets]))
                display(self.optimizer.bo_progress)
                 
            with self.outputopt_panel:
                print('Output')
                # display(widgets.VBox([self.inputopt_panel, self.optimizer.runcalc_button]))

            display(widgets.VBox([self.inputopt_panel, self.accordion_outopt]))


   
    def display_import_panel(self):

        with self.import_panel:

            upload_text = "Upload Structure file"
            uploadWidget = widgets.HTML(value=f"<b><font color='black'>{upload_text}</b>")

            selectBox = widgets.HBox([self.__selectcif_button, self.__uploadcif_button])
            uploadBox = widgets.VBox([uploadWidget, selectBox])
            display(uploadBox)
            

    def display_config_panel(self):

        """Display the building panel for Builder Module"""

        with self.config_panel:

            text_calculator = 'Calculator: '
            calcWidget = widgets.HTML(value=f"<b><font color='black'>{text_calculator}</b>")

            calculator = widgets.HBox([calcWidget, self.calc_buttons])
            # display(configWidget)

            text_remote = 'Remote Cluster Configuration:'
            remoteWidget = widgets.HTML(value=f"<b><font color='black'>{text_remote}</b>")

            with self.remote_panel:
                text_cluster = ' If necessary, add your own remote cluster: '
                textWidget = widgets.HTML(value=f"<b><font color='black'>{text_cluster}</b>")

                display(widgets.VBox([widgets.HBox([self.remote.clusterslist_dropdown, self.remote.password_text,
                                                    self.remote.passphrase_text, self.remote.cluster_button]),
                                      textWidget,
                                      widgets.HBox(
                                          [self.remote.remote_text, self.remote.host_text, self.remote.user_text,
                                           self.remote.port_text]),
                                      widgets.HBox([self.remote.id_checkbox, self.remote.upload_button]),
                                      self.remote.remote_button]))

            text_vasp = 'VASP Configuration:'
            vaspWidget = widgets.HTML(value=f"<b><font color='black'>{text_vasp}</b>")

            with self.vaspconfig_panel:

                display(widgets.VBox([self.vaspconfig.accuracy_buttons,
                                      self.vaspconfig.xc_buttons,
                                      self.vaspconfig.setup_buttons,
                                      self.vaspconfig.settings_text]))

            

            text_job = 'Job file:'
            jobWidget = widgets.HTML(value=f"<b><font color='black'>{text_job}</b>")

            #with self.jobconfig_panel:
            self.accordion_jobtext = widgets.Accordion(children=[self.vaspconfig.jobfile_text], selected_index=None)
            self.accordion_jobtext.set_title(0, 'VASP job file')

            text_qe = 'Quantum Espresso Configuration:'
            qeWidget = widgets.HTML(value=f"<b><font color='black'>{text_qe}</b>")

            self.accordion_qejobtext = widgets.Accordion(children=[self.qeconfig.jobfile_text], selected_index=None)
            self.accordion_qejobtext.set_title(0, 'Quantum Espresso job file')

                #display(widgets.VBox([self.accordion_jobtext]))

            display(widgets.VBox([calculator,
                                  remoteWidget,
                                  self.remote_panel,
                                  vaspWidget,
                                  self.vaspconfig_panel,
                                  #jobWidget,
                                  self.accordion_jobtext,
                                  qeWidget,
                                  self.qeconfig_panel,
                                  self.accordion_qejobtext]))

    def display_jobs_panel(self):
        with self.jobs_panel:
            show_buttons = widgets.HBox([self.remote.boupdate_button, self.remote.jobsrunning_button,
                                         self.remote.bojob_text,self.remote.outputjob_button,
                                         self.remote.restartjob_button, self.remote.killjob_button])
            display(widgets.VBox([show_buttons, self.jobslist_panel]))

    def display_data_panel(self):
        with self.data_panel:
            show_buttons = widgets.HBox([self.database.showdb_button, self.database.db_buttons])
            display(widgets.VBox([show_buttons, self.showdb_panel]))

    def __create_display_button(self):
        """Build the Display button widget."""
        
        # create widget
        self.__display_button = LoadedButton(description='Display', value=None, disabled=True)
        
        # update function
        self.__display_button.on_click(self.__display_button_on_clicked)

    def __display_button_on_clicked(self, b):
        
        # self.search_menu.phasediagram_button.value = None
        
        b.value = [str(self.search_menu.search_select.value.split()[0]), self.search_menu.viewer_radio.value,
                   self.search_menu.unitcell_radio.value]
        
    def pd_clicked(self, b):

        self.__display_button.value = ['Phase Diagram',  self.search_menu.structure.value]
        
        with self.error_panel:
            clear_output()
            print('Button Clicked !')
            print(self.__display_button.value)

    def update_slab_button(self, *args):

        if self.orientation_ckecked:
            self.builder.slab_button.disabled = False
            self.builder.slab_button.style.button_color = 'lightgreen'
        else:
            self.builder.slab_button.disabled = True
            self.builder.slab_button.button_style = 'danger'

    def update_orientation_button(self, *args):

        if self.search_menu.molecule_checkbox.value:
            self.builder.orientation_button.disabled = True
            self.builder.orientation_button.button_style = 'danger'
        else:
            if str(self.search_menu.search_select.value.split()[0]) == 'None':
                self.builder.orientation_button.disabled = True
                self.builder.orientation_button.button_style = 'danger'

            else:
                self.builder.orientation_button.disabled = False
                self.builder.orientation_button.style.button_color = 'lightgreen'

    def jobsrunning_clicked(self, b):
        with self.jobslist_panel:
            clear_output()
            self.remote.runningjobs_print()

    def boupdate_clicked(self, b):

        def takeSecond(elem):
            return int(elem[4:])

        boids = {}
        with open('../asapdata/boids.json', mode='r', encoding='utf-8') as fd:
            json_file = json.load(fd)
            boids = json_file

        list_bo = sorted(list(boids.keys()), key=takeSecond, reverse=True)
        with self.jobslist_panel:
            clear_output()

            for boid in list_bo:

                if boids[boid]['status'] != 'killed':

                    if not os.path.isfile('../asapdata/{0}/output.txt'.format(boid)):

                        text_bo = boid + ' (pid: {0})'.format(boids[boid]['pid'])
                        bo_label = widgets.HTML(value=f"<font color='black'>{text_bo}")
                        text_system = boids[boid]['system']
                        system_label = widgets.HTML(value=f"<b><font color='black'>{text_system}</b>")

                        bo_progress = widgets.IntProgress(value=0, min=0, max=15, step=1,
                                                          description='progress: ')

                        energy_label = widgets.HTML(value=f"<b><font color='black'>{' Minimum Energy:'}</b>")
                        value_label = widgets.HTML(value=f"<b><font color='black'>{''}</b>")
                        bojob_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                                      height='50px', overflow_y='auto'))
                        display(bojob_panel)
                        with bojob_panel:
                            display(widgets.HBox([system_label, bo_label, bo_progress, energy_label, value_label]))

                    else:

                        filebo = open('../asapdata/{0}/output.txt'.format(boid), 'r')
                        lines = filebo.readlines()

                        if len(lines) >= 1:

                            if lines[-1] == 'Finished':

                                min_energy = 0
                                for line in lines[:-1]:
                                    energy = float(line.split()[3])
                                    if energy < min_energy:
                                        min_energy = energy
                                        list_optinput = line.split()[5:]

                                last_step = lines[-2]
                                steps = last_step.split()[1]
                                current_step = int(steps.split('/')[0]) + 1
                                max_step = int(steps.split('/')[1])
                                init_step = last_step.split()[2]
                                init_step = init_step.split(')')[0]
                                init_step = int(init_step[1:])

                                text_bo = boid + ' (pid: {0})'.format(boids[boid]['pid'])
                                bo_label = widgets.HTML(value=f"<font color='black'>{text_bo}")
                                text_system = boids[boid]['system']
                                system_label = widgets.HTML(value=f"<b><font color='black'>{text_system}</b>")

                                total_steps = max_step + init_step
                                bo_progress = widgets.IntProgress(value=current_step, min=0, max=total_steps, step=1,
                                                                  description='progress: ', bar_style='success')

                                energy_label = widgets.HTML(value=f"<b><font color='black'>{' Minimum Energy:'}</b>")
                                value_label = widgets.HTML(value=f"<font color='black'>{str(min_energy)+ ' eV'}")
                                vector = ''
                                for optvalue in list_optinput:
                                    vector += str(optvalue) + ' '

                                input_label = widgets.HTML(value=f"<b><font color='black'>{'Best parameters: '}</b>")
                                vector_label = widgets.HTML(value=f"<font color='black'>{vector}")

                                bojob_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                                              height='90px', overflow_y='auto'))
                                display(bojob_panel)
                                with bojob_panel:
                                    display(widgets.HBox([system_label, bo_label, bo_progress, energy_label, value_label]))
                                    display(widgets.HBox([input_label, vector_label]))

                            else:

                                min_energy = 0
                                for line in lines:
                                    energy = float(line.split()[3])
                                    if energy < min_energy:
                                        min_energy = energy

                                last_step = lines[-1]
                                steps = last_step.split()[1]
                                current_step = int(steps.split('/')[0]) + 1
                                max_step = int(steps.split('/')[1])
                                init_step = last_step.split()[2]
                                init_step = init_step.split(')')[0]
                                init_step = int(init_step[1:])

                                text_bo = boid + ' (pid: {0})'.format(boids[boid]['pid'])
                                bo_label = widgets.HTML(value=f"<font color='black'>{text_bo}")
                                text_system = boids[boid]['system']
                                system_label = widgets.HTML(value=f"<b><font color='black'>{text_system}</b>")

                                total_steps = max_step+init_step
                                bo_progress = widgets.IntProgress(value=current_step, min=0, max=total_steps, step=1,
                                                                  description='progress: ')

                                energy_label = widgets.HTML(value=f"<b><font color='black'>{' Minimum Energy:'}</b>")
                                value_label = widgets.HTML(value=f"<b><font color='black'>{str(min_energy) + ' eV'}</b>")
                                bojob_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                                                   height='50px', overflow_y='auto'))
                                display(bojob_panel)
                                with bojob_panel:
                                    display(widgets.HBox([system_label, bo_label, bo_progress, energy_label, value_label]))

                        else:

                            text_bo = boid + ' (pid: {0})'.format(boids[boid]['pid'])
                            bo_label = widgets.HTML(value=f"<font color='black'>{text_bo}")
                            text_system = boids[boid]['system']
                            system_label = widgets.HTML(value=f"<b><font color='black'>{text_system}</b>")

                            bo_progress = widgets.IntProgress(value=0, min=0, max=15, step=1,
                                                              description='progress: ')

                            energy_label = widgets.HTML(value=f"<b><font color='black'>{' Minimum Energy:'}</b>")
                            value_label = widgets.HTML(value=f"<b><font color='black'>{''}</b>")
                            bojob_panel = ipywidgets.Output(layout=Layout(width='99.6%', border='solid 1px',
                                                                          height='50px', overflow_y='auto'))
                            display(bojob_panel)
                            with bojob_panel:
                                display(widgets.HBox([system_label, bo_label, bo_progress, energy_label, value_label]))

            # self.remote.boupdate_print()

    def booutput_clicked(self, b):
        with self.ssh_panel:
            clear_output()
            jobid = self.remote.bojob_text.value.strip()

            image_convergence = '../asapdata/{0}/bo_convergence.png'.format(jobid)

            from IPython.display import Image, display

            display(Image(filename=image_convergence))


            image_acquisition = '../asapdata/{0}/bo_acquisition.png'.format(jobid)

            if os.path.isfile(image_acquisition):
                display(Image(filename=image_acquisition))



    def killjob_clicked(self, b):
        jobid = self.remote.bojob_text.value.strip()

        boids = {}
        with open('../asapdata/boids.json', mode='r', encoding='utf-8') as fd:
            json_file = json.load(fd)
            boids = json_file

            boids[jobid]['status'] = 'killed'

        with open('../asapdata/boids.json', mode='w', encoding='utf-8') as fp:
            json.dump(boids, fp)

            os.system('kill -9 ' + str(boids[jobid]['pid']))

        with self.ssh_panel:
            clear_output()
            print('BO {0} killed !'.format(jobid))

    def showdb_clicked(self, b):
        with self.showdb_panel:
            clear_output()
            self.database.pprint()

    def checkorien_clicked(self, b):

        self.all_slabs, _, self.all_orientations, _ = \
            self.builder.check_orientation(str(self.search_menu.search_select.value.split()[0]))
        b.value = [self.all_orientations, self.all_slabs]

        self.builder.orientation_select.options = self.all_orientations
        self.builder.orientation_select.value = str(self.all_orientations[0])

        # self.builder.slabs_dropdown.options = self.all_orientations
        # self.builder.slabs_dropdown.value = self.all_orientations[0]

        self.orientation_ckecked = True

        self.builder.slab_button.disabled = False
        self.builder.slab_button.button_style = 'success'

        with self.error_panel:
            clear_output()
            print('Orientation Clicked !')
            print(self.all_orientations)
            print(len(self.all_slabs))
            # print(b.value[1])
            # print(b.value[2])

    def nanoparticle_clicked(self, b):

        mp_id = str(self.search_menu.search_select.value.split()[0])

        self.np_atoms = self.nanoparticle.create_nanoparticle(mp_id)

        b.value = [self.np_atoms]

        self.__display_button.value = ['Nanoparticle',  self.np_atoms]

        with self.error_panel:
            clear_output()
            print('Nanoparticle Clicked !')
            print(mp_id)
            print(self.np_atoms)
            print(self.nanoparticle.list_energies)
            print(self.nanoparticle.list_families)
            print(self.__display_button.value)


    def rundft_clicked(self, b):
        print('Bonjour')

        self.aimd.rundft_action(self.__display_button.value, self.remote.password, self.remote.passphrase)

    def runoptmd_clicked(self, b):
        print('Bonjour')

        self.optmd.runoptmd_action(self.__display_button.value, self.remote.password, self.remote.passphrase)

    def slab_clicked(self, b):

        id_slab = int(self.builder.orientation_select.value.split('-')[0])

        with self.error_panel:
            clear_output()
            print('Slab Clicked !')
            print(self.all_slabs[id_slab])
            print(self.builder.orientation_select.value)

        self.slab_chosen = self.builder.create_slab(str(self.search_menu.search_select.value.split()[0]),
                                                    self.builder.orientation_select.value,
                                                    self.all_slabs[id_slab])

        b.value = [self.slab_chosen]

        self.__display_button.value = ['Slab',  
                                       self.search_menu.structure.value, 
                                       self.builder.miller_text.value,
                                       self.builder.orientation_select.value,
                                       self.builder.slabsize_text.value,
                                       self.builder.vacuum_text.value,
                                       self.builder.slabcenter_checkbox.value,
                                       self.builder.supercellx_text.value,
                                       self.builder.supercelly_text.value,
                                       self.builder.supercellz_text.value
                                      ]

        # self.builder.slab_button.disabled = True
        # self.builder.slab_button.button_style = 'warning'

        with self.sites_panel:
            clear_output()
            structure = self.builder.slab_button.value[0]
            self.builder.adsorbate_sites(self.slab_chosen)

            print("Relevant adsorption sites: ", len(self.builder.ads_sites['all']))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax = plot_slab(structure, ax, adsorption_sites=True)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.show()

        with self.error_panel:
            # clear_output()
            print('Slab Clicked !')
            # print(orientation_selected)
            # print(self.builder.orientation_button.value[2])
            # print(self.builder.orientation_button.value[0][orientation_selected])
            # print(self.slab_chosen)
            # print(b.value[0])
            
            print('{0} slabs for {1}'.format(len(self.builder.slab_button.value),
                                             str(self.search_menu.search_select.value.split()[0])))

            slab = AseAtomsAdaptor.get_atoms(self.builder.slab_button.value[0])

            print(slab.get_cell())
            # print(self.builder.slab_button.value[0])
            
            # print(self.builder.ads_sites['all'][0])
            # print(self.builder.ads_sites['all'][1])

    def checkmol_clicked(self, b):
        
        from ase.collections import g2
        
        try:
            dropdown = []
            
            if self.builder.adatom_text.value in extra:
                dropdown.append('{0} in ase database'.format(self.builder.adatom_text.value))
                
            # results = self.search_menu.search_molecules(self.builder.adatom_text.value)
            # b.value = results
            elif self.builder.adatom_text.value in g2.names:
                dropdown.append('{0} in ase database'.format(self.builder.adatom_text.value))

            # for molecule in results:
                # dropdown.append(molecule['task_id'] + ' {1} ({0:.4f})'.format(molecule['IE'], molecule['formula']))
                
        except:
            dropdown = []
            self.builder.adsorbate_button.disabled = True
            self.builder.adsorbate_button.button_style = 'warning'

        with self.adsorb_panel:
            
            text_cell = 'Results for {0}:'.format(self.builder.adatom_text.value)
            self.molText.value = f"<b><font color='black'>{text_cell}</b>"
            self.builder.molad_dropdown.disabled = False
            
            if len(dropdown) == 0:
                self.builder.molad_dropdown.options = [('No molecule found', 0)]
                self.builder.adsorbate_button.disabled = True
                self.builder.adsorbate_button.style.button_color = 'red'
            else:
                self.builder.molad_dropdown.options = dropdown
                self.builder.molad_dropdown.value = dropdown[0]
                self.builder.adsorbate_button.disabled = False
                self.builder.adsorbate_button.style.button_color = 'lightgreen'

    def update_adsorbate_button(self, *args):
        
        if self.builder.adatom_radio.value == 'Atom':
            self.builder.adsorbate_button.disabled = False

    def bo_calculate(self, x):
        
        # print(x[0])
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
        if self.builder.adatom_radio.value == 'Atom':
            
            structure = self.builder.add_adsorbate_ase(x=x_input, y=y_input, h=h_input)

        elif self.builder.adatom_radio.value == 'Molecule':
            from ase.build import molecule
            atoms = molecule(self.builder.adatom_text.value)

            structure = self.builder.add_adsorbate_mol(atoms, x=x_input, y=y_input, h=h_input,
                                                       alpha=alpha_input, beta=beta_input, gamma=gamma_input)

        # current_structure = AseAtomsAdaptor.get_atoms(structure)

        # self.optimizer.runcalc_button.value = current_structure

        # print('VASP Calculator / adsorbate')
        # print(self.current_structure)
        # print(current_structure)
        # print(current_structure.symbols)
        # bo_step = self.optimizer.bo_progress.value
        # self.createinput(atoms=current_structure, bo_step=bo_step)
        # print('Remote workdir: ', self.remote_workdir)
        # self.remote.execute_job(path=self.remote_workdir, system=str(current_structure.symbols), structure='bo_'+str(bo_step))

        # No need to change the button value ? Send directly the structure to optimizer function
        # self.optimizer.calculate(method=self.calc_buttons.value)
        # y_next = self.optimizer.e_slab
        # print(x, y_next)
        
        # self.optimizer.bo_progress.value += 1

        structure = AseAtomsAdaptor.get_atoms(structure)
        self.optimizer.runcalc_button.value = structure
        # No need to change the button value ? Send directly the structure to optimizer function
        self.optimizer.calculate()
        y_next = self.optimizer.e_slab
        print(x, y_next)

        self.optimizer.bo_progress.value += 1

        return y_next

    def createinput_clicked(self, b, path_input='input_vasp/'):

        with self.outputopt_panel:
            clear_output()

            if self.calc_buttons.value == 'EMT':
                print('Input files for EMT no necesary')
            elif self.calc_buttons.value == 'VASP':
                print(self.__display_button.value)
                if self.__display_button.value[0] == 'Adsorbate':
                    print(self.__display_button.value[9]+'_on_'+self.__display_button.value[1]
                          + '_' + self.__display_button.value[2].replace(" ", "")
                          + '_' + self.__display_button.value[3]
                          + '_' + self.__display_button.value[4])
                elif self.__display_button.value[0] == 'Slab':
                    print(self.__display_button.value[1]
                          + '_' + self.__display_button.value[2].replace(" ", "")
                          + '_' + self.__display_button.value[3]
                          + '_' + self.__display_button.value[4])
                elif self.__display_button.value[0] == 'Optimize':
                    print(self.__display_button.value[1]['symbols'])

                self.vaspconfig.create_input(prec=self.vaspconfig.accuracy_buttons.value)
                print('Input files for VASP created. Available in input_vasp repertory')
                import os
                self.remote.create_workdir(path='ASAP/testASAP2')
                self.remote.put_files(path='ASAP/testASAP2', path_input=path_input)
                self.remote.command_ls(path='ASAP/testASAP2')

    def createinput(self, atoms=None, path_input='input_vasp/', bo_step=None, bulk=False, ssh=True):
        with self.outputopt_panel:
            # clear_output()
            if bulk:  # Bulk input

                if self.calc_buttons.value == 'EMT':
                    print('Input files for EMT no necesary')
                elif self.calc_buttons.value == 'VASP':
                    workdir = str(atoms.symbols)

                    self.vaspconfig.create_input(atoms=atoms, prec=self.vaspconfig.accuracy_buttons.value)

                    print('Input files for VASP created. Available in input_vasp repertory')
                    import os

                    self.remote_workdir = 'ASAP/Bulk/' + workdir

                    self.remote.create_workdir(path=self.remote_workdir)
                    self.remote.put_files(path=self.remote_workdir, path_input=path_input)

            else:  # Slab input
                if self.calc_buttons.value == 'EMT':
                    print('Input files for EMT no necesary')
                elif self.calc_buttons.value == 'VASP':
                    print(self.__display_button.value)
                    if self.__display_button.value[0] == 'Adsorbate':
                        workdir = self.__display_button.value[9]+'_on_'+self.__display_button.value[1]\
                                  + '_' + self.__display_button.value[2].replace(" ", "")\
                                  + '_' + self.__display_button.value[3]\
                                  + '_' + self.__display_button.value[4]
                        print(workdir)

                        if bo_step is None:
                            self.remote_workdir = 'ASAP/Adsorbate/' + workdir
                        else:
                            self.remote_workdir = 'ASAP/Adsorbate/' + workdir + '/' + bo_step

                    elif self.__display_button.value[0] == 'Slab':

                        workdir = self.__display_button.value[1]\
                                  + '_' + self.builder.orientation_select.value.replace(" ", "")\
                                  + '_' + self.__display_button.value[4]\
                                  + '_' + self.__display_button.value[5]
                        print(workdir)

                        if bo_step is None:
                            self.remote_workdir = 'ASAP/Slab/' + workdir
                        else:
                            self.remote_workdir = 'ASAP/Slab/' + workdir + '/' + bo_step

                    elif self.__display_button.value[0] == 'Optimize':
                        print(self.__display_button.value[1]['symbols'])

                    self.vaspconfig.create_input(atoms=atoms, prec=self.vaspconfig.accuracy_buttons.value)
                    print('Input files for VASP created. Available in input_vasp repertory')
                    import os
                    if ssh:
                        self.remote.create_workdir(path=self.remote_workdir)
                        self.remote.put_files(path=self.remote_workdir, path_input=path_input)
                    # self.remote.command_ls(path=self.remote_workdir)

    def bulk_calc(self, mpid=None):

        if mpid not in self.bulk_db.keys():
            self.vaspconfig.bulk_input(mpid=mpid)
            bulk_structure = self.rester.get_structure_by_material_id(mpid)
            structure = self.rester.get_entry_by_material_id(mpid)

            bulk_atoms = AseAtomsAdaptor.get_atoms(bulk_structure)

            # self.optimizer.runcalc_button.value = bulk_atoms
            self.builder.sernergy_button.value = bulk_atoms

            self.createinput(atoms=bulk_atoms, bo_step=None, bulk=True)
            print('Remote workdir: ', self.remote_workdir)
            system = str(bulk_atoms.symbols)
            self.remote.execute_job(path=self.remote_workdir, system=system, structure='bulk')

            self.bulk_db['bulk'][mpid] = structure.energy

            # self.database.add_entry(mpid=mpid, energy=structure.energy, structure=None, bulk=True, slab=False)
            self.database.write_to_json()

        else:
            print('Bulk Energy already in database for ', mpid)
            return self.bulk_db['bulk'][mpid]

    def runcalc(self):
        with self.outputopt_panel:
            clear_output()
            if self.calc_buttons.value == 'EMT':

                if self.__display_button.value[0] == 'Adsorbate' or self.__display_button.value[0] == 'Optimize':
                    # print(self.builder.ads_slab)
                    structure = AseAtomsAdaptor.get_structure(self.builder.ads_slab)
                    self.current_structure = AseAtomsAdaptor.get_atoms(structure)
                    self.optimizer.runcalc_button.value = self.current_structure
                    self.optimizer.calculate('onecalc')

                    y_next = self.optimizer.e_slab

                    print('Energy: ', y_next)
                else:
                    print('No Run calculation for single slab with EMT')

            elif self.calc_buttons.value == 'VASP':
                if self.__display_button.value[0] == 'Slab':
                    print('VASP Calculator / slab')
                    # print(self.current_structure[0])
                    self.current_structure = AseAtomsAdaptor.get_atoms(self.current_structure[0])
                    print(self.current_structure)
                    print(self.current_structure.symbols)
                    self.optimizer.runcalc_button.value = self.current_structure

                    self.createinput(atoms=self.current_structure, bo_step=None)

                    print('Remote workdir: ', self.remote_workdir)
                    system = str(self.current_structure.symbols)
                    self.remote.execute_job(path=self.remote_workdir, system=system, structure='slab')

                elif self.__display_button.value[0] == 'Adsorbate':
                    print('VASP Calculator / adsorbate')
                    # print(self.current_structure)
                    self.current_structure = AseAtomsAdaptor.get_atoms(self.current_structure)
                    print(self.current_structure)
                    print(self.current_structure.symbols)

                    self.createinput(atoms=self.current_structure, bo_step=None)
                    print('Remote workdir: ', self.remote_workdir)
                    system = str(self.current_structure.symbols)
                    self.remote.execute_job(path=self.remote_workdir, system=system, structure='adsorbate')

                elif self.__display_button.value[0] == 'Optimize':

                    print('VASP Calculator / After Optimize')
                    # print(self.current_structure)

                    print(self.current_structure)
                    print(self.current_structure.symbols)

    def surfacenergy_clicked(self, b):

        if self.calc_buttons.value == 'EMT':

            with self.ssh_panel:
                clear_output()
                print('EMT')
                print(type(self.slab_chosen))
                print(self.slab_chosen)

        elif self.calc_buttons.value == 'VASP':

            with self.ssh_panel:

                print('VASP Calculator / slab')

                # Check if Bulk calculation done. If not, do the bulk calculation
                mpid = str(self.search_menu.search_select.value.split()[0])
                self.bulk_calc(mpid=mpid)
                current_structure = AseAtomsAdaptor.get_atoms(self.slab_chosen)
                print(current_structure.symbols)
                print('Dropdown: ', str(self.search_menu.search_select.value.split()[0]))

                # self.optimizer.runcalc_button.value = current_structure
                self.builder.sernergy_button.value = current_structure
                self.createinput(atoms=current_structure, bo_step=None)
                print('Remote workdir: ', self.remote_workdir)
                self.remote.execute_job(path=self.remote_workdir, system=str(current_structure.symbols),
                                        structure='slab')

    def runcalc_clicked(self, b):

        with self.outputopt_panel:
            clear_output()
            if self.calc_buttons.value == 'EMT':

                if self.__display_button.value[0] == 'Adsorbate' or self.__display_button.value[0] == 'Optimize':
                    # print(self.builder.ads_slab)
                    structure = AseAtomsAdaptor.get_structure(self.builder.ads_slab)
                    current_structure = AseAtomsAdaptor.get_atoms(structure)
                    self.optimizer.runcalc_button.value = current_structure
                    self.optimizer.calculate('onecalc')

                    y_next = self.optimizer.e_slab

                    print('Energy: ', y_next)
                else:
                    print('No Run calculation for single slab with EMT')

            elif self.calc_buttons.value == 'VASP':

                if self.__display_button.value[0] == 'Slab':
                    print('VASP Calculator / slab')

                    # Check if Bulk calculation done. If not, do the bulk calculation
                    mpid = str(self.search_menu.search_select.value.split()[0])
                    self.bulk_calc(mpid=mpid)

                    # print(self.current_structure[0])

                    current_structure = AseAtomsAdaptor.get_atoms(self.current_structure[0])
                    print(current_structure)
                    print(current_structure.symbols)

                    print('Dropdown: ', str(self.search_menu.search_select.value.split()[0]))

                    self.optimizer.runcalc_button.value = current_structure
                    self.createinput(atoms=current_structure, bo_step=None)
                    print('Dropdown: ', str(self.search_menu.search_select.value.split()[0]))
                    print('Remote workdir: ', self.remote_workdir)
                    self.remote.execute_job(path=self.remote_workdir, system=str(current_structure.symbols))

                elif self.__display_button.value[0] == 'Adsorbate':

                    import subprocess
                    print('Run outside process')
                    print(os.getcwd())
                    cmd = 'python src/bayesopt.py'
                    p = subprocess.Popen(cmd)
                    print('Outside process running')

                    print('VASP Calculator / adsorbate')
                    # print(self.current_structure)
                    current_structure = AseAtomsAdaptor.get_atoms(self.current_structure)
                    print(current_structure)
                    print(current_structure.symbols)

                    self.createinput(atoms=current_structure, bo_step=None)
                    print('Remote workdir: ', self.remote_workdir)
                    self.remote.execute_job(path=self.remote_workdir, system=str(current_structure.symbols),
                                            structure='surface')

                elif self.__display_button.value[0] == 'Optimize':

                    print('VASP Calculator / After Optimize')
                    # print(self.current_structure)

                    print(self.current_structure)
                    print(self.current_structure.symbols)

                # self.optimizer.calculate('onecalc')

    def runopt_vasp(self):

        # Create a new id for the bayesian optimization
        # All the bo id are saved inside a json file.

        from ase import Atoms

        domain_dico = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                       'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                       'h': {'name': 'h', 'type': 'continuous', 'domain': (1, 5)},
                       'alpha': {'name': 'alpha', 'type': 'continuous', 'domain': (-180, 180)},
                       'beta': {'name': 'beta', 'type': 'continuous', 'domain': (-180, 180)},
                       'gamma': {'name': 'gamma', 'type': 'continuous', 'domain': (-180, 180)}}

        domain_dico_discrete = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                                'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                                'h': {'name': 'h', 'type': 'continuous', 'domain': (1, 5)},
                                'alpha': {'name': 'alpha', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'beta': {'name': 'beta', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'gamma': {'name': 'gamma', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)}}

        self.emin = None
        self.slabmin = None
        self.optimizer.bo_progress.value = 0
        self.optimizer.bo_progress.max = self.optimizer.maxbo_slider.value + self.optimizer.ninit

        slab = AseAtomsAdaptor.get_atoms(self.builder.slab_button.value[0])

        xads = str(self.builder.ads_sites['all'][self.builder.sites_dropdown.value][0])
        yads = str(self.builder.ads_sites['all'][self.builder.sites_dropdown.value][1])

        workdirsystem = self.__display_button.value[9] + '_on_' + self.__display_button.value[1] \
                  + '_' + self.__display_button.value[2].replace(" ", "") \
                  + '_' + self.__display_button.value[3] \
                  + '_' + self.__display_button.value[4]

        boids = {}

        if os.path.isfile('../asapdata/boids.json'):

            with open('../asapdata/boids.json', mode='r', encoding='utf-8') as fd:
                json_file = json.load(fd)
                boids = json_file

                max_key = 0
                for key in list(boids.keys()):
                    if int(key[4:]) > max_key:
                        max_key = int(key[4:])

                boid = max_key + 1
                workdir = 'boid' + str(boid)
                boids[workdir] = {'system': workdirsystem, 'pid': None, 'status': 'running', 'step': -1}

            with open('../asapdata/boids.json', mode='w', encoding='utf-8') as fp:
                json.dump(boids, fp)

        else:
            boid = 0
            workdir = 'boid' + str(boid)
            boids[workdir] = {'system': workdirsystem, 'pid': None, 'status': 'running', 'step': -1}

            with open('../asapdata/boids.json', mode='w', encoding='utf-8') as fp:
                json.dump(boids, fp)

        local_workdir = '../asapdata/' + workdir
        os.makedirs(local_workdir + '/vasp/', exist_ok=True)

        if self.calc_buttons.value == 'VASP':
            # Run outside Bayesian Optimization
            with self.outputopt_panel:
                clear_output()
                print(self.__display_button.value[0])
                print('Test BO to run a job on cluster with VASP')

                with open(local_workdir+'/domain', 'w') as file:

                    file.write(self.builder.adatom_radio.value)
                    file.write('\n')
                    for parameter in self.optimizer.parameters_select.value:
                        file.write(str(parameter))
                        file.write('\n')

                ase.io.write(local_workdir+'/slab', slab, format='vasp')

                if self.builder.adatom_radio.value == 'Atom':

                    molecule = self.builder.adatom_text.value
                    molecule = Atoms(molecule, positions=[(0., 0., 0.)])

                    ase.io.write(local_workdir+'/molecule', molecule, format='xyz')

                elif self.builder.adatom_radio.value == 'Molecule':

                    from ase.build import molecule
                    atoms = molecule(self.builder.adatom_text.value)

                    ase.io.write(local_workdir+'/molecule', atoms, format='xyz')

                current_structure = AseAtomsAdaptor.get_atoms(self.current_structure)
                print(current_structure)
                os.makedirs(local_workdir+'/vasp/', exist_ok=True)

                self.createinput(atoms=current_structure,  path_input=local_workdir+'/vasp/', bo_step=None, ssh=False)
                import shutil
                shutil.copy2('input_vasp/POTCAR', local_workdir+'/vasp')  # complete target filename given
                shutil.copy2('input_vasp/KPOINTS', local_workdir+'/vasp')  # complete target filename given
                shutil.copy2('input_vasp/INCAR', local_workdir+'/vasp')  # complete target filename given
                shutil.copy2('input_vasp/job_submit', local_workdir+'/vasp')  # complete target filename given
                shutil.copy2('input_vasp/submit_file.sh', local_workdir+'/vasp')  # complete target filename given

                import subprocess
                print('Run outside process')
                print(os.getcwd())
                cmd = 'python src/bayesopt.py'

                maxbo = str(self.optimizer.maxbo_slider.value)
                initbo = str(self.optimizer.initbo_slider.value)


                if self.remote.password is None and self.remote.passphrase is not None:
                    cmd = 'python src/bayesopt.py' + ' --passphrase ' + self.remote.passphrase \
                          + ' --workdir ' + workdir + ' --maxbo ' + maxbo + ' --initbo ' + initbo + ' --x ' + xads \
                          + ' --y ' + yads
                elif self.remote.password is not None and self.remote.passphrase is None:
                    cmd = 'python src/bayesopt.py' + ' --password ' + self.remote.password + ' --workdir ' + workdir \
                          + ' --maxbo ' + maxbo + ' --initbo ' + initbo + ' --x ' + xads + ' --y ' + yads
                elif self.remote.password is not None and self.remote.passphrase is not None:
                    cmd = 'python src/bayesopt.py' + ' --password ' + self.remote.password + \
                           ' --passphrase ' + self.remote.passphrase + ' --workdir ' + workdir + ' --maxbo ' + maxbo \
                          + ' --initbo ' + initbo + ' --x ' + xads + ' --y ' + yads

                print(cmd)
                p = subprocess.Popen(cmd, shell=True)
                # print(p.pid)

                boids = {}

                with open('../asapdata/boids.json', mode='r', encoding='utf-8') as fd:
                    json_file = json.load(fd)
                    boids = json_file

                    boids[workdir]['pid'] = p.pid

                with open('../asapdata/boids.json', mode='w', encoding='utf-8') as fp:
                    json.dump(boids, fp)

                print('Bayesian optimization running with id {0}'.format(workdir))

    def runopt_clicked(self, b):

        domain_dico = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                       'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                       'h': {'name': 'h', 'type': 'continuous', 'domain': (1, 5)},
                       'alpha': {'name': 'alpha', 'type': 'continuous', 'domain': (-180, 180)},
                       'beta': {'name': 'beta', 'type': 'continuous', 'domain': (-180, 180)},
                       'gamma': {'name': 'gamma', 'type': 'continuous', 'domain': (-180, 180)}}

        domain_dico_discrete = {'x': {'name': 'x', 'type': 'continuous', 'domain': (0, 0.99)},
                                'y': {'name': 'y', 'type': 'continuous', 'domain': (0, 0.99)},
                                'h': {'name': 'h', 'type': 'continuous', 'domain': (1, 5)},
                                'alpha': {'name': 'alpha', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'beta': {'name': 'beta', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)},
                                'gamma': {'name': 'gamma', 'type': 'discrete', 'domain': np.arange(-180, 190, 30)}}

        self.emin = None
        self.slabmin = None
        self.optimizer.bo_progress.value = 0
        self.optimizer.bo_progress.max = self.optimizer.maxbo_slider.value + self.optimizer.ninit

        slab = AseAtomsAdaptor.get_atoms(self.builder.slab_button.value[0])

        if self.calc_buttons.value == 'EMT':

            with self.outputopt_panel:
                clear_output()

                self.domain = []
                # print("HELLO")

                for parameter in self.optimizer.parameters_select.value:
                    # print(parameter)
                    # self.domain.append(domain_dico[parameter])
                    self.domain.append(domain_dico_discrete[parameter])

                # print(self.domain)
                # --- Solve your problem

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
                                              initial_design_numdata=self.optimizer.ninit,
                                              acquisition_type='EI',
                                              exact_feval=True)

                myBopt.run_optimization(max_iter=self.optimizer.maxbo_slider.value)

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

                # for i in range(10):

                # structure = self.builder.add_adsorbate_ase(h=h_lists[i])
                # b.value = structure
                # print('h = {0}'.format(h_lists[i]))
                # self.optimizer.calculate()
                # if self.emin == None or self.emin > self.optimizer.e_slab:
                # self.emin = self.optimizer.e_slab
                # atoms = AseAtomsAdaptor.get_atoms(structure)
                # self.slabmin = atoms
            with self.outputopt2_panel:
                clear_output()
                myBopt.plot_acquisition()
                myBopt.plot_convergence()
                # print(np.round(myBopt.x_opt, 2))
                print('Optimum parameters: ', np.round(myBopt.x_opt, 2))
                # print('Minimum Energy: ', np.round(myBopt.y_opt, 2))

            if self.builder.adatom_radio.value == 'Atom':

                structure = self.builder.add_adsorbate_ase(x=x_input, y=y_input, h=h_input)
            elif self.builder.adatom_radio.value == 'Molecule':
                from ase.build import molecule
                atoms = molecule(self.builder.adatom_text.value)

                structure = self.builder.add_adsorbate_mol(atoms, x=x_input, y=y_input, h=h_input, alpha=alpha_input,
                                                           beta=beta_input, gamma=gamma_input)

            atoms = AseAtomsAdaptor.get_atoms(structure)
            self.slabmin = atoms
            self.current_structure = atoms

            self.__display_button.value = ['Optimize',
                                           self.slabmin]

        elif self.calc_buttons.value == 'VASP':

            self.runopt_vasp()

    def add_adsorbate_clicked(self, b):
        
        slab = AseAtomsAdaptor.get_atoms(self.builder.slab_button.value[0])
            
        if self.builder.adatom_radio.value == 'Atom':
            # Adsorbate is an atom
            
            self.optimizer.parameters_select.options = ['x', 'y', 'h']
            self.optimizer.parameters_select.values = ['x', 'y', 'h']
            
            b.value = self.builder.add_adsorbate_ase()
            self.current_structure = b.value

        else:
            # Adsorbate is a molecule
            
            self.optimizer.parameters_select.options = ['x', 'y', 'h', 'alpha', 'beta', 'gamma']
            self.optimizer.parameters_select.values = ['x', 'y', 'h', 'alpha', 'beta', 'gamma']

            from ase.build import molecule
            atoms = molecule(self.builder.adatom_text.value)

            with self.error_panel:
                clear_output()
                print('Adsorbate Clicked !')

                print('{0} adsorbate for {1}'.format(atoms,
                                                     str(self.search_menu.search_select.value.split()[0])))
           
            # print(AseAtomsAdaptor.get_atoms(molecule))
            b.value = self.builder.add_adsorbate_mol(atoms)
            self.current_structure = b.value

        self.__display_button.value = ['Adsorbate',
                                       self.search_menu.structure.value,
                                       self.builder.miller_text.value,
                                       self.builder.slabsize_text.value,
                                       self.builder.vacuum_text.value,
                                       self.builder.slabcenter_checkbox.value,
                                       self.builder.supercellx_text.value,
                                       self.builder.supercelly_text.value,
                                       self.builder.supercellz_text.value,
                                       self.builder.adatom_text.value,
                                       self.builder.sites_dropdown.value,
                                       self.builder.molad_dropdown.value,
                                       self.builder.height_slider.value,
                                       self.builder.alpha_slider.value,
                                       self.builder.beta_slider.value,
                                       self.builder.gamma_slider.value]
             
        with self.error_panel:
            clear_output()
            print('Adsorbate Added !')
            
            print('{0} adsorbate for {1}'.format(self.builder.adatom_text.value,
                                                 str(self.search_menu.search_select.value.split()[0])))

    # Function to update the state of Display button
    # This should be activatedd/deactivated everytime the value of Search_Button becomes False/None/something
    
    def update_display_button(self, *args):
        
        if self.search_menu.search_button.value is not False and not None:
            self.__display_button.button_style = 'warning'
            self.__display_button.tooltip = 'Select the output to display'
            
            if self.search_menu.viewer_radio.value is not None:
                self.__display_button.tooltip = 'Ready to display'
                self.__display_button.disabled = False
                self.__display_button.style.button_color = 'lightgreen'
        
    # Function to update the state of Display button
    # This should be activatedd/deactivated everytime the value of Search_Button becomes False/None/something
    
    def update_description_panel(self, *args):
        
        if str(self.search_menu.search_select.value.split()[0]) is None:

            # Description_Label.value = 'Waiting for structure'
            with self.description_panel:
                clear_output()
                print('Waiting for structure')

        else:

            if self.search_menu.molecule_checkbox.value:

                molecule = self.search_menu.get_molecule(str(self.search_menu.search_select.value.split()[0]))
                with self.description_panel:
                    clear_output()
                    print('Molecule Properties:')
                    print('Formula: ', molecule['formula'])
                    print('Ionization Energy: ', molecule['IE'])
                    if 'EA'  in molecule.keys() :
                        print('Affinity Energy: ', molecule['EA'])

            else:

                structure = self.rester.get_structure_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                conventional_structure = \
                    SpacegroupAnalyzer(structure).get_conventional_standard_structure(international_monoclinic=True)
                
                with self.description_panel:
                    clear_output()
                    # print(str(self.search_menu.search_select.value.split()[0]))
                    if self.search_menu.unitcell_radio.value == 'primitive':
                        print(structure)
                    else:
                        print(conventional_structure)
    
    # Function to update the state of Output panel
    # This should be activate everytime the value of Search Select or Unit Cell or Display Output Select changes                    
    def update_output_panel(self, *args):
        
        with self.output_panel:
            clear_output()
            
            if self.__display_button.value[0] == 'Phase Diagram':
                self.__savecif_button.disabled = True
                self.__savecif_button.button_style = 'warning'
                self.__saveposcar_button.disabled = True
                self.__saveposcar_button.button_style = 'warning'

                # print('Phase Diagram !!')
                pd = PhaseDiagram(self.search_menu.search_button.value[1])
                plotterpd = PDPlotter(pd)   
                plotterpd.show()   
            
            elif self.__display_button.value[0] == 'Input file':
                # structure_ase = atoms
                #atoms = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
                atoms = read(self.__display_button.value[1])
                try:
                    display(viewer_mod.view_ngl(atoms))
                except:
                    display(view(structure_ase, viewer='x3d'))


            elif self.__display_button.value[1] == '3D structure': 

                if self.search_menu.unitcell_radio.value == 'primitive':
                    # print(Search_Button.value)
                    structure = self.rester.get_structure_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                    atoms = AseAtomsAdaptor.get_atoms(structure)
                    entry_mp = self.rester.get_entry_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                    structure_ase = atoms
                    
                    try:
                        display(viewer_mod.view_ngl(atoms))
                    except:
                        display(view(structure_ase, viewer='x3d'))
                        
                    # viewstruct = nglview.show_pymatgen(structure)
                    # viewase = nglview.show_ase(atoms)
                    # display(viewase)
                    # display(viewstruct)
                    # display(view(atoms, viewer=viewer))

                elif self.search_menu.unitcell_radio.value == 'conventional':
                    structure = self.rester.get_structure_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                    conventional_structure = \
                        SpacegroupAnalyzer(structure).get_conventional_standard_structure(international_monoclinic=True)
                    atoms = AseAtomsAdaptor.get_atoms(conventional_structure)
                    entry_mp = self.rester.get_entry_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                    structure_ase = atoms
                    
                    try:
                        display(viewer_mod.view_ngl(atoms))
                    except:
                        display(view(structure_ase, viewer='x3d'))
                    # viewstruct = nglview.show_pymatgen(conventional_structure)
                    # viewase = nglview.show_ase(atoms)
                    # display(viewase)

                    # display(viewstruct)
                    # display(view(atoms, viewer=viewer))
                
                self.__savecif_button.value = atoms
                self.__savecif_button.disabled=False
                self.__savecif_button.button_style = 'success'
                self.__saveposcar_button.disabled = False
                self.__saveposcar_button.button_style = 'success'

                
            elif self.__display_button.value[1] == 'Bands structure':
                self.__savecif_button.disabled = True
                self.__savecif_button.button_style = 'warning'
                self.__saveposcar_button.disabled = True
                self.__saveposcar_button.button_style = 'warning'
                from pymatgen.electronic_structure.plotter import BSPlotter
                bs = self.rester.get_bandstructure_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                plotter = BSPlotter(bs)
                plotter.show()

            elif self.__display_button.value[1] == 'DOS':
                self.__savecif_button.disabled = True
                self.__savecif_button.button_style = 'warning'
                self.__saveposcar_button.disabled = True
                self.__saveposcar_button.button_style = 'warning'

                from pymatgen.electronic_structure.plotter import DosPlotter
                dos = self.rester.get_dos_by_material_id(str(self.search_menu.search_select.value.split()[0]))
                dos_plotter = DosPlotter()
                dos_plotter.add_dos_dict(dos.get_spd_dos())
                dos_plotter.show()

            elif self.__display_button.value[0] == 'Nanoparticle':

                atoms = self.np_atoms

                try:
                    display(viewer_mod.view_ngl(atoms))
                except:
                    display(view(atoms, viewer='x3d'))


                self.__savecif_button.value = atoms
                self.__savecif_button.disabled = False
                self.__savecif_button.button_style = 'success'
                self.__saveposcar_button.disabled = False
                self.__saveposcar_button.button_style = 'success'


            elif self.__display_button.value[0] == 'Slab':
                
                atoms = AseAtomsAdaptor.get_atoms(self.slab_chosen)
                
                try:
                    display(viewer_mod.view_ngl(atoms))
                except:
                    display(view(atoms, viewer='x3d'))
                        
                self.__savecif_button.value = atoms
                self.__savecif_button.disabled = False
                self.__savecif_button.button_style = 'success'
                self.__saveposcar_button.disabled = False
                self.__saveposcar_button.button_style = 'success'

            elif self.__display_button.value[0] == 'Adsorbate':

                # structure = self.builder.adsorbate_button.value[self.builder.sites_dropdown.value]
                # atoms = AseAtomsAdaptor.get_atoms(structure)
                # display(viewer_mod.view_ngl(atoms))
                # print(self.builder.ads_slab)
                
                try:
                    display(viewer_mod.view_ngl(self.builder.ads_slab))
                except:
                    display(view(self.builder.ads_slab, viewer='x3d'))
                
                self.__savecif_button.value = self.builder.ads_slab
                
                self.__savecif_button.disabled = False
                self.__savecif_button.button_style = 'success'
                self.__saveposcar_button.disabled = False
                self.__saveposcar_button.button_style = 'success'
                
            elif self.__display_button.value[0] == 'Optimize':

                try:
                    display(viewer_mod.view_ngl(self.slabmin))
                except:
                    display(view(self.slabmin, viewer='x3d'))
                
                self.__savecif_button.value = self.slabmin
                
                self.__savecif_button.disabled = False
                self.__savecif_button.button_style = 'success'
                self.__saveposcar_button.disabled = False
                self.__saveposcar_button.button_style = 'success'

        # with self.adsorb_panel:

            # if self.__display_button.value[0] == 'Adsorbate':
                # clear_output()
                # print("Relevant adsorption sites: ")
     
                # fig = plt.figure()
                # ax = fig.add_subplot(111)
                # ax = plot_slab(structure, ax, adsorption_sites=True)
                # ax.set_xticks([])
                # ax.set_yticks([])
                # plt.show()
                # plot_slab(structure, ax, adsorption_sites=False, decay=0.09)

    def update_sites_panel(self, *args):
        
        with self.sites_panel:

            if self.__display_button.value[0] == 'Slab':
                clear_output()
                
                structure = self.builder.slab_button.value[0] 
                
                self.builder.adsorbate_sites(self.current_structure)
                # print(n_sites)
                
                print("Relevant adsorption sites: ", len(self.builder.ads_sites['all']))
     
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax = plot_slab(structure, ax, adsorption_sites=True)
                ax.set_xticks([])
                ax.set_yticks([])
                plt.show()
                    

extra = {
    'Be2': {
        'symbols': 'BeBe',
        'positions': [[0, 0, 1.0106],
                      [0, 0, -1.0106]]},
    'C7NH5': {
        'symbols': 'C7NH5',
        'positions': [[-1.593581, -1.142601, 0.],
                      [-2.235542, 0.095555, 0.],
                      [-0.204885, -1.210726, 0.],
                      [0.549645, -0.025355, 0.],
                      [1.976332, -0.085321, 0.],
                      [-0.099258, 1.220706, 0.],
                      [-1.488628, 1.273345, 0.],
                      [3.136871, -0.128138, 0.],
                      [-2.177996, -2.060896, 0.],
                      [-3.323594, 0.141242, 0.],
                      [0.301694, -2.173705, 0.],
                      [0.488716, 2.136782, 0.],
                      [-1.987765, 2.240495, 0.]]},
    'BDA': {
        # 1,4-Benzodiamine
        # aka p-Aminoaniline; p-Benzenediamine; p-Diaminobenzene;
        #     p-Phenylenediamine; Paraphenylen-diamine
        # PBE-gpaw relaxed
        'symbols': 'C6H4N2H4',
        'positions': [[0.004212, 1.406347, 0.061073],
                      [1.193490, 0.687096, 0.029481],
                      [1.190824, -0.690400, -0.028344],
                      [0.000295, -1.406191, -0.059503],
                      [-1.186974, -0.685668, -0.045413],
                      [-1.185376, 0.690203, 0.009452],
                      [2.147124, 1.219997, 0.064477],
                      [2.141593, -1.227477, -0.054266],
                      [-2.138408, -1.222814, -0.095050],
                      [-2.137740, 1.226930, 0.023036],
                      [-0.006314, 2.776024, 0.186278],
                      [-0.007340, -2.777839, -0.159936],
                      [0.844710, -3.256543, 0.110098],
                      [-0.854965, -3.253324, 0.130125],
                      [0.845826, 3.267270, -0.055549],
                      [-0.854666, 3.254654, -0.092676]]},
    'biphenyl': {
        # PBE-gpaw relaxed
        'symbols': 'C6H5C6H5',
        'positions': [[-0.74081, -0.00000, -0.00003],
                      [-1.46261, -1.20370, -0.00993],
                      [-2.85531, -1.20350, -0.00663],
                      [-3.55761, -0.00000, -0.00003],
                      [-2.85531, 1.20350, 0.00667],
                      [-1.46261, 1.20370, 0.00997],
                      [-0.92071, -2.14850, 0.00967],
                      [-3.38981, -2.15110, -0.00083],
                      [-4.64571, -0.00000, -0.00003],
                      [-3.38981, 2.15110, 0.00077],
                      [-0.92071, 2.14850, -0.00963],
                      [3.55849, -0.00000, -0.00003],
                      [2.85509, -0.86640, -0.83553],
                      [1.46289, -0.87000, -0.83153],
                      [0.73969, -0.00000, -0.00003],
                      [1.46289, 0.87000, 0.83157],
                      [2.85509, 0.86640, 0.83547],
                      [4.64659, -0.00000, -0.00003],
                      [3.39189, -1.53770, -1.50253],
                      [0.91869, -1.53310, -1.50263],
                      [0.91869, 1.53310, 1.50267],
                      [3.39189, 1.53770, 1.50257]]},
    'C60': {
        # Buckminsterfullerene, I*h symm.
        # The Buckyball has two degrees of freedom, the C-C bond, and the
        # C=C bond. This is an LDA-gpaw relaxed structure with bond lengths
        # 1.437 and 1.385.
        # Experimentally, the two bond lengths are 1.45 and 1.40 Angstrom.
        'symbols': 'C60',
        'positions': [[2.2101953, 0.5866631, 2.6669504],
                      [3.1076393, 0.1577008, 1.6300286],
                      [1.3284430, -0.3158939, 3.2363232],
                      [3.0908709, -1.1585005, 1.2014240],
                      [3.1879245, -1.4574599, -0.1997005],
                      [3.2214623, 1.2230966, 0.6739440],
                      [3.3161210, 0.9351586, -0.6765151],
                      [3.2984981, -0.4301142, -1.1204138],
                      [-0.4480842, 1.3591484, 3.2081020],
                      [0.4672056, 2.2949830, 2.6175264],
                      [-0.0256575, 0.0764219, 3.5086259],
                      [1.7727917, 1.9176584, 2.3529691],
                      [2.3954623, 2.3095689, 1.1189539],
                      [-0.2610195, 3.0820935, 1.6623117],
                      [0.3407726, 3.4592388, 0.4745968],
                      [1.6951171, 3.0692446, 0.1976623],
                      [-2.1258394, -0.8458853, 2.6700963],
                      [-2.5620990, 0.4855202, 2.3531715],
                      [-0.8781521, -1.0461985, 3.2367302],
                      [-1.7415096, 1.5679963, 2.6197333],
                      [-1.6262468, 2.6357030, 1.6641811],
                      [-3.2984810, 0.4301871, 1.1204208],
                      [-3.1879469, 1.4573895, 0.1996030],
                      [-2.3360261, 2.5813627, 0.4760912],
                      [-0.5005210, -2.9797771, 1.7940308],
                      [-1.7944338, -2.7729087, 1.2047891],
                      [-0.0514245, -2.1328841, 2.7938830],
                      [-2.5891471, -1.7225828, 1.6329715],
                      [-3.3160705, -0.9350636, 0.6765268],
                      [-1.6951919, -3.0692581, -0.1976564],
                      [-2.3954901, -2.3096853, -1.1189862],
                      [-3.2214182, -1.2231835, -0.6739581],
                      [2.1758234, -2.0946263, 1.7922529],
                      [1.7118619, -2.9749681, 0.7557198],
                      [1.3130656, -1.6829416, 2.7943892],
                      [0.3959024, -3.4051395, 0.7557638],
                      [-0.3408219, -3.4591883, -0.4745610],
                      [2.3360057, -2.5814499, -0.4761050],
                      [1.6263757, -2.6357349, -1.6642309],
                      [0.2611352, -3.0821271, -1.6622618],
                      [-2.2100844, -0.5868636, -2.6670300],
                      [-1.7726970, -1.9178969, -2.3530466],
                      [-0.4670723, -2.2950509, -2.6175105],
                      [-1.3283500, 0.3157683, -3.2362375],
                      [-2.1759882, 2.0945383, -1.7923294],
                      [-3.0909663, 1.1583472, -1.2015749],
                      [-3.1076090, -0.1578453, -1.6301627],
                      [-1.3131365, 1.6828292, -2.7943639],
                      [0.5003224, 2.9799637, -1.7940203],
                      [-0.3961148, 3.4052817, -0.7557272],
                      [-1.7120629, 2.9749122, -0.7557988],
                      [0.0512824, 2.1329478, -2.7937450],
                      [2.1258630, 0.8460809, -2.6700534],
                      [2.5891853, 1.7227742, -1.6329562],
                      [1.7943010, 2.7730684, -1.2048262],
                      [0.8781323, 1.0463514, -3.2365313],
                      [0.4482452, -1.3591061, -3.2080510],
                      [1.7416948, -1.5679557, -2.6197714],
                      [2.5621724, -0.4853529, -2.3532026],
                      [0.0257904, -0.0763567, -3.5084446]]}}                
