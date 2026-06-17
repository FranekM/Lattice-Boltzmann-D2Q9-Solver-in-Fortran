import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import f90nml

fps = 24

config_names = [
    'swd_amp_plots',
    'dwd_amp_plots',
    'swd_vis_plots',
    'swd_prf_plots',
    'dwd_prf_plots',
    'swd_fld_anim',
    'dwd_fld_anim',
    'swd_prf_anim',
    'swd_prf_anim'
]
configs_raw = f90nml.read("../../config.in")['experiment']
config = {}
for c in configs_raw:
    c = c['cfg']
    config[c['name']] = c

def frame_fld_reader(cfg):
    reader = pd.read_csv(
            '../../' + cfg["filename"],
            sep       = '\t',
            names     = ['t', 'y', 'x', 'rho', 'uy', 'ux'],
            chunksize = cfg["Nx"] * cfg["Ny"]
    )

    for frame in reader:
        yield frame

def frame_prf_reader(cfg):
    frames = pd.read_csv(
            '../../' + cfg["filename"],
            sep         = '\t',
            names       = ['t', 'x', 'rho', 'uy', 'ux'],
            chunksize   = cfg["Nx"]
    )
    for frame in frames:
        yield frame


def anim_prf(cfg):
    frames = frame_prf_reader(cfg)
    first = next(frames)

    fig, ax = plt.subplots()
    line, = ax.plot(first[cfg["disp_x"]], first[cfg["disp_y"]])

    def update(frame):
        line.set_data(frame[cfg["disp_x"]], frame[cfg["disp_y"]])
        return line,

    a = anim.FuncAnimation(
            fig, update,
            frames   = frames,
            interval = 1000.0 / fps,
            blit     = False
    )
    plt.show()

def anim_fld(cfg):
    frames = frame_fld_reader(cfg)
    first = next(frames)

    fig, ax = plt.subplots()

    U = first['ux'].values.reshape(cfg["Ny"], cfg["Nx"]).T
    V = first['uy'].values.reshape(cfg["Ny"], cfg["Nx"]).T
    R = first['rho'].values.reshape(cfg["Ny"], cfg["Nx"]).T

    x, y = np.arange(cfg["Nx"]), np.arange(cfg["Ny"])
    X, Y = np.meshgrid(x, y)
    qv = ax.quiver(
            X, Y, U, V,
            scale       = 0.001,
            scale_units = 'xy',
            linewidth   = 2.0,
            color       = 'black'
    )

    im = ax.imshow(R, cmap='Blues_r')
    fig.colorbar(im, ax = ax)
    ax.invert_yaxis()

    def update(frame):
        U = frame["ux"].values.reshape(cfg["Ny"], cfg["Nx"]).T
        V = frame["uy"].values.reshape(cfg["Ny"], cfg["Nx"]).T
        R = frame["rho"].values.reshape(cfg["Ny"], cfg["Nx"]).T

        qv.set_UVC(U, V)
        im.set_data(R)
        return qv, im

    a = anim.FuncAnimation(
            fig, update,
            frames   = frames,
            interval = 1000.0 / fps,
            blit     = False
    )
    plt.show()


def plot_prf(cfg):
    df = pd.read_csv(
            '../../' + cfg["filename"],
            sep         = '\t',
            names       = ['t', 'x', 'rho', 'uy', 'ux']
    )

    for t, group in df.groupby(['t']):
        plt.plot(group[cfg["disp_x"]], group[cfg["disp_y"]], label=f"time step: {t}")
        plt.xlabel('x')
        plt.ylabel('value')

    plt.legend()
    plt.show()

def plot_amp(cfg, omega_idx=None):
    df = pd.read_csv(
            '../../' + cfg["filename"],
            sep         = '\t',
            names       = ['omega', 't', 'A']
    )

    if omega_idx is not None:
        values = sorted(df.omega.unique())
        selected = [values[i] for i in omega_idx]
        df = df[df.omega.isin(selected)]

    for omega, group in df.groupby(['omega']):
        plt.plot(group[cfg["disp_x"]], group[cfg["disp_y"]], label=f"omega: {omega}")
        plt.xlabel('time step')
        plt.ylabel('amplitude')

    plt.legend()
    plt.show()


def plot_vis(cfg):
    df = pd.read_csv(
            '../../' + cfg["filename"],
            sep='\t',
            names=['omega', 't', 'A']
    )

    omegas = []
    measured_nus = []
    expected_nus = []

    for omega, group in df.groupby("omega"):
        t = group.t
        A = group.A
        slope, offset = np.polyfit(t, np.log(A), 1)

        Ny = 100
        k = 2 * np.pi / Ny
        nu_measured = -slope / k ** 2
        nu_expected = (1.0 / omega - 0.5) / 3.0

        omegas.append(omega)
        measured_nus.append(nu_measured)
        expected_nus.append(nu_expected)

    plt.plot(omegas, measured_nus, label="measured nu")
    plt.plot(omegas, expected_nus, label="expected nu")
    plt.xlabel('omega')
    plt.ylabel('nu')
    plt.legend()
    plt.grid()
    plt.show()


anim_prf(config['dwd_prf_anim'])
anim_prf(config['swd_prf_anim'])
anim_fld(config['dwd_fld_anim'])
anim_fld(config['swd_fld_anim'])
plot_prf(config['dwd_prf_plots'])
plot_prf(config['swd_prf_plots'])
plot_amp(config['dwd_amp_plots'])
plot_amp(config['swd_amp_plots'])
plot_vis(config['swd_vis_plots'])