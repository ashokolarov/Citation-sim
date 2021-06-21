from DataProcessing import DynamicData
import numpy as np
from constants import *
from matplotlib import pyplot as plt
from model import Model


def sym_mode_comparison(label):
    if label == 'sp':
        t_num, dat_num = model.response_symmetric_mode('sp', 10)
        sym_mode_fl = flight.SymMode(1, 1, 0)
        fl = sym_mode_fl[0]
        title = 'Short period'
    else:
        t_num, dat_num = model.response_symmetric_mode('ph', 120)
        sym_mode_fl = flight.SymMode(1, 1, 0)
        fl = sym_mode_fl[1]
        title = 'Phugoid'

    fig, axs = plt.subplots(2, 2, sharex='col')

    axs[0, 0].ticklabel_format(useOffset=False)
    axs[0, 0].ticklabel_format(style='plain')

    axs[0, 0].plot(t_num, dat_num[0], label='Numerical data', color=num_color, linewidth=2, linestyle=num_style)
    axs[0, 0].plot(fl[0], fl[1], label='Flight data', color=flight_color, linewidth=2, linestyle=flight_style)
    axs[0, 0].set_ylabel('Velocity[m/s]', fontsize=18)
    axs[0, 0].grid(True)

    axs[0, 1].plot(t_num, dat_num[1], color=num_color, linewidth=2, linestyle=num_style)
    axs[0, 1].plot(fl[0], fl[2], color=flight_color, linewidth=2, linestyle=flight_style)
    axs[0, 1].set_ylabel(r'Angle of attack $\alpha$[deg]', fontsize=18)
    axs[0, 1].grid(True)

    axs[1, 0].plot(t_num, dat_num[2], color=num_color, linewidth=2, linestyle=num_style)
    axs[1, 0].plot(fl[0], fl[3], color=flight_color, linewidth=2, linestyle=flight_style)
    axs[1, 0].set_ylabel(r'Pitch $\theta$[deg]', fontsize=18)
    axs[1, 0].set_xlabel('Time[s]', fontsize=18)
    axs[1, 0].grid(True)

    axs[1, 1].plot(t_num, dat_num[3], color=num_color, linewidth=2, linestyle=num_style)
    axs[1, 1].plot(fl[0], fl[4],  color=flight_color, linewidth=2, linestyle=flight_style)
    axs[1, 1].set_ylabel('Pitch rate q[deg/s]', fontsize=18)
    axs[1, 1].set_xlabel('Time[s]', fontsize=18)
    axs[1, 1].grid(True)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    fig.suptitle(f'Comparison between model and flight data during {title} motion', size='x-large', weight='bold')
    fig.legend(lines, labels, loc='lower center', fancybox=True, shadow=True, ncol=2)
    plt.show()


def asym_mode_comparison(label):
    if label == 'dr':
        t_num, dat_num = model.response_asymmetric_mode('dr', 15)
        asym_mode_fl = flight.AsymMode(1, 1, 1, 1)
        fl = asym_mode_fl[0]
        title = 'Dutch roll'

    elif label == 'apr':
        t_num, dat_num = model.response_asymmetric_mode('apr', 15)
        asym_mode_fl = flight.AsymMode(1, 1, 1, 1)
        fl = asym_mode_fl[2]
        title = 'Aperiodic roll'

    else:
        t_num, dat_num = model.response_asymmetric_mode('spi', 50)
        asym_mode_fl = flight.AsymMode(1, 1, 1, 1)
        fl = asym_mode_fl[3]
        title = 'Spiral'

    fig, axs = plt.subplots(2, 2, sharex='col')

    axs[0, 0].plot(t_num, dat_num[0], label='Numerical data', color=num_color, linewidth=2, linestyle=num_style)
    axs[0, 0].plot(fl[0], fl[2], label='Flight data', color=flight_color, linewidth=2, linestyle=flight_style)
    axs[0, 0].set_ylabel(fr'Sideslip $\beta$[deg]', fontsize=18)
    axs[0, 0].grid(True)

    axs[0, 1].plot(t_num, dat_num[3], color=num_color, linewidth=2, linestyle=num_style)
    axs[0, 1].plot(fl[0], fl[5], color=flight_color, linewidth=2, linestyle=flight_style)
    axs[0, 1].set_ylabel('Yaw rate r[deg/s]', fontsize=18)
    axs[0, 1].grid(True)

    axs[1, 0].plot(t_num, dat_num[1], color=num_color, linewidth=2, linestyle=num_style)
    axs[1, 0].plot(fl[0], fl[3], color=flight_color, linewidth=2, linestyle=flight_style)
    axs[1, 0].set_ylabel('Roll angle φ[deg]', fontsize=18)
    axs[1, 0].set_xlabel('Time[s]', fontsize=18)
    axs[1, 0].grid(True)

    axs[1, 1].plot(t_num, dat_num[2], color=num_color, linewidth=2, linestyle=num_style)
    axs[1, 1].plot(fl[0], fl[4], color=flight_color, linewidth=2, linestyle=flight_style)
    axs[1, 1].set_ylabel('Roll rate p[deg/s]', fontsize=18)
    axs[1, 1].set_xlabel('Time[s]', fontsize=18)
    axs[1, 1].grid(True)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

    fig.suptitle(f'Comparison between model and flight data during {title} motion', size='x-large', weight='bold')
    fig.legend(lines, labels, loc='lower center', fancybox=True, shadow=True, ncol=2)
    plt.show()


def snakiness():
    t, dat_num = model.response_asymmetric_mode('dr', 15)
    asym_mode_fl = flight.AsymMode(1, 1, 1, 1)
    fl = asym_mode_fl[0]
    title = 'Snakiness Plot'

    fig, axs = plt.subplots(1, 1)

    axs.plot(dat_num[3], dat_num[2], label='Numerical', color=num_color, linewidth=2, linestyle=num_style)
    axs.plot(fl[5], fl[4], label='Flight data', color=flight_color, linewidth=2, linestyle=flight_style)
    axs.set_ylabel('Roll rate p[deg/s]', fontsize=18)
    axs.set_xlabel('Yaw rate r[deg/s]', fontsize=18)
    axs.set_aspect('equal', 'box')
    axs.grid(True)

    fig.tight_layout()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    model = Model()
    flight = DynamicData()

    sym_mode_comparison('sp')
    sym_mode_comparison('ph')


