from pymatgen import Structure, Lattice, MPRester, Molecule
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from matplotlib import pyplot as plt


from fireworks import Workflow

from atomate.vasp.fireworks.core import StaticFW



# Note that you must provide your own API Key, which can
# be accessed via the Dashboard at materialsproject.org
mpr = MPRester('4oBTKz0pkFSg9EUQ')

from fireworks import LaunchPad
lpad = LaunchPad.auto_load()
#lpad.reset('', require_password=False)

from atomate.vasp.workflows.presets.core import wf_static
from atomate.vasp.workflows.base.adsorption import get_wf_slab
from atomate.vasp.powerups import add_additional_fields_to_taskdocs, add_tags
struct = mpr.get_structure_by_material_id("mp-23") # fcc Ni
struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
slabs = generate_all_slabs(struct, 1, 5.0, 2.0, center_slab=True)
slab_dict = {slab.miller_index:slab for slab in slabs}

ni_slab_111 = slab_dict[(1, 1, 1)]
print(type(ni_slab_111))
print(ni_slab_111)

db_file='/home/f-gimbert/atomate/config/db.json'

fw = StaticFW(structure=ni_slab_111, db_file=db_file)

wf = Workflow([fw],name=' Slab static wf from BOCALC')
#wf = wf_static(ni_slab_111)
#wf = get_wf_slab(ni_slab_111, db_file='/home/f-gimbert/atomate/config/db.json')
#new_wf = add_additional_fields_to_taskdocs(wf, {"system":"surface"})
#new_wf = add_tags(wf, ['surface'])
lpad.add_wf(wf)