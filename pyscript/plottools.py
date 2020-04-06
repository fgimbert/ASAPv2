#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math


def plot_esm1(filename, titlename=''):
    esm1 = np.loadtxt(filename, comments='#')

    fig = plt.figure(
        figsize=(6, 4),  # inch
        dpi=100,  # dpi
        edgecolor='black',
        linewidth='1'
    )

    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.suptitle(titlename)

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_xlabel('z (A)')
    ax1.set_ylabel('rho (e/A)')
    ax2.set_xlabel('z (A)')
    ax2.set_ylabel('V_hartree (eV)')
    ax3.set_xlabel('z (A)')
    ax3.set_ylabel('V_ion (eV)')
    ax4.set_xlabel('z (A)')
    ax4.set_ylabel('V_electrostatic (eV)')

    ax4.axhline(0.0, linewidth='1', linestyle='dashed', color='black')

    ax1.plot(esm1[:, 0], esm1[:, 1], color='black', linestyle='solid')
    ax2.plot(esm1[:, 0], esm1[:, 2], color='black', linestyle='solid')
    ax3.plot(esm1[:, 0], esm1[:, 3], color='black', linestyle='solid')
    ax4.plot(esm1[:, 0], esm1[:, 4], color='black', linestyle='solid')

    plt.show()


def plot_1drism(filename, titlename='', normalization=False, max_x=None):
    with open(filename, 'r') as file:
        data_1drism = file.readlines()
        # line_molecules = data_1drism[3].split()
        line_atoms = data_1drism[4].split()

    data_1drism = np.loadtxt(filename, comments='#', skiprows=5)

    number_of_sublines = int(math.ceil((len(line_atoms) - 1) / 3))
    number_of_subplots = len(line_atoms) - 1

    fig = plt.figure(
        figsize=(15, 2 * number_of_sublines),  # inch
        dpi=100,  # dpi
        edgecolor='black',
        linewidth='1'
    )

    fig.subplots_adjust(wspace=0.5, hspace=0.75)
    fig.suptitle(titlename)

    factor_norm = 1

    for n_plot in range(1, number_of_subplots + 1):
        if normalization:
            factor_norm = data_1drism[-1, n_plot]

        ax1 = fig.add_subplot(number_of_sublines, 3, n_plot)
        if max_x is not None:
            ax1.set_xlim([0, max_x])

        ax1.set_xlabel('r (A)')
        ax1.set_ylabel('rdf')

        if (number_of_subplots % 3) == 0:
            if n_plot <= (number_of_subplots // 3 - 1) * 3:
                ax1.set_xlabel('')
            # ax1.set_xticklabels([])

        ax1.set_title(line_atoms[n_plot])
        ax1.plot(data_1drism[:, 0], data_1drism[:, n_plot] / factor_norm, color='black', linestyle='solid')

    #    plt.savefig('%s.pdf'%titlename)
    plt.show()


def plot_rism1(filename, titlename=''):
    with open(filename, 'r') as file:
        data_rism1 = file.readlines()
        line_data = data_rism1[1].split()

    data_rism1 = np.loadtxt(filename, comments='#', skiprows=2)

    line_data = line_data[1:]
    new_line_data = [line_data[0] + ' ' + line_data[1]]
    for i in range(2, len(line_data), 3):
        new_line_data.append(line_data[i] + ' ' + line_data[i + 1] + ' ' + line_data[i + 2])

    number_of_sublines = int(math.ceil((len(new_line_data) - 1) / 3))
    number_of_subplots = len(new_line_data) - 1

    fig = plt.figure(
        figsize=(15, 2 * number_of_sublines),  # inch
        dpi=100,  # dpi
        edgecolor='black',
        linewidth='1'
    )

    fig.subplots_adjust(wspace=0.5, hspace=0.75)
    fig.suptitle(titlename)

    for n_plot in range(1, number_of_subplots + 1):
        ax1 = fig.add_subplot(number_of_sublines, 3, n_plot)
        ax1.set_xlabel(new_line_data[0])
        ax1.set_ylabel(new_line_data[n_plot])
        # ax1.set_title(line_atoms[n_plot])
        ax1.plot(data_rism1[:, 0], data_rism1[:, n_plot], color='black', linestyle='solid')

    plt.show()
