#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pyscript.localfiles import LocalFiles
from ase.io import read
from ase.io.espresso import read_espresso_extra
import os


project = 'Steve'
prj_info = {
    'project': project,
    'homedir': os.path.expanduser('~')+'/Program',
    'prefixes':  ['CuCOmu'],
}

rjob = LocalFiles(prj_info)

scf_result = {}
extra_result = {}
for prefix in rjob.prefixes:
    scf_result[prefix] = read(rjob.srcout[prefix], index=':', format='espresso-out')
    with open(rjob.srcout[prefix], mode='rt') as fd:
        extra_result[prefix] = read_espresso_extra(fd, index=':')

for prefix in rjob.prefixes:
    extra = extra_result[prefix]
    print('***** Total charge from input file *****')
    print("tot_charge = ", f"{extra[-1]['tot_charge']:>12.5}")
    if extra[-1]['trism'] is True:
        print('***** Chemical potential of solvation from 1D-RISM *****')
        print('                          Closure      Gauss Fluctuation')
        for solvx in extra[-1]['mols']:
            for idy, solvy in enumerate(extra[-1]['mols']):
                if idy == 0:
                    print(f"{solvx:^10}{solvy:^10}{extra[-1]['1d_chem_pot'][solvx][solvy][0]:>16.6f}"
                          f"{extra[-1]['1d_chem_pot'][solvx][solvy][1]:>16.6f}")
                else:
                    print(f"          {solvy:^10}{extra[-1]['1d_chem_pot'][solvx][solvy][0]:>16.6f}"
                          f"{extra[-1]['1d_chem_pot'][solvx][solvy][1]:>16.6f}")
        print('***** Solvent information from 3D/Laue-RISM *****')
        print("Image  solvent    # of solv.    Solv. chg(e).   Chem. pot. (kcal/mol)")
        for idx, image in enumerate(extra[0:-1]):
            for idy, solv in enumerate(extra[-1]['mols']):
                print(f"{idx+1:^5}{solv:^10}{image['solvent_num'][idy]:>12.6f}{image['solvent_chg'][idy]:>16.6f}"
                      f"{image['solvent_pot'][idy]:18.6f}")
    if extra[-1]['lfcp'] is True:
        print('***** FCP information *****')
        for image in extra[0:-1]:
            print("Total charge", f"{image['fcp_current_chg']:12.5f}", " --(next)-->", f"{image['fcp_next_chg']:>12.5f}")
            print("Fermi energy", f"{image['fcp_fermi']:12.5f}", "(eV)")
            print("Force on FCP", f"{image['fcp_force']:12.5f}", "(eV)")
