from matplotlib import rc
rc("animation", html="jshtml", embed_limit=200)
from IPython.display import Video, display

from matplotlib.animation import FFMpegWriter
import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import f90nml
import math
import re

from pathlib import Path


fps = 24


def init(config_file):
    path = Path(f"../config/{config_file}.in")
    # print(path.absolute())
    configs_raw = f90nml.read(path)['experiment']
    config = {}
    if len(configs_raw) > 1:
        for c in configs_raw:
            c = c['cfg']
            config[c['name']] = c
    else:
        c = configs_raw['cfg']
        config[c['name']] = c
    return config


def _frame_generator(filename):
    """
    Zwraca kolejne DataFrame'y odpowiadające jednemu timestepowi.
    """

    reader = pd.read_csv(
        filename,
        sep = "\t",
        names = ["omega", "t", "y", "x", "rho", "uy", "ux"],
        chunksize = 4096,
    )

    buffer = pd.DataFrame()
    for chunk in reader:
        buffer = pd.concat([buffer, chunk], ignore_index=True)
        while len(buffer):
            t = buffer.iloc[0]["t"]
            mask = buffer["t"] == t
            # cały bufor należy do jednego timestepu
            if mask.all():
                break
            yield buffer[mask].copy()
            buffer = buffer[~mask].reset_index(drop=True)

    if len(buffer):
        yield buffer

def log_frame_reader(cfg, home_dir):
    if home_dir is None:
        home_dir = "../"
    path = Path(home_dir) / cfg["log"]["file"]
    ext = cfg["log"]["ext"]
    stem = path.name

    pattern = re.compile( rf"{re.escape(stem)}(\d+){re.escape(ext)}$" )
    readers = []

    for file in sorted(path.parent.iterdir()):
        if not pattern.match(file.name):
            continue
        if file.stat().st_size == 0:
            continue
        readers.append(_frame_generator(file))

    if not readers:
        raise FileNotFoundError(f"No log files found for {stem}")

    while True:
        frames = []
        for reader in readers:
            try:
                frames.append(next(reader))
            except StopIteration:
                pass
        if not frames:
            break
        yield (
            pd.concat(frames, ignore_index=True)
              .sort_values(["y", "x"])
              .reset_index(drop=True)
        )

def log_dataframe_reader(cfg, home_dir):
    if home_dir is None:
        home_dir = "../"
    path = Path(home_dir) / cfg["log"]["file"]
    ext = cfg["log"]["ext"]
    stem = path.name

    pattern = re.compile( rf"{re.escape(stem)}(\d+){re.escape(ext)}$" )
    frames = []

    for file in sorted(path.parent.iterdir()):
        if not pattern.match(file.name):
            continue
        if file.stat().st_size == 0:
            continue
        frames.append(
            pd.read_csv(
                file,
                sep = "\t",
                names = ["omega", "t", "y", "x", "rho", "uy", "ux"],
            )
        )

    if not frames:
        raise FileNotFoundError(f"No log files found for {stem}")

    return (
        pd.concat(frames, ignore_index=True)
          .sort_values(["omega", "t", "y", "x"])
          .reset_index(drop=True)
    )


def field_animation(cfg, qv_step, qv_scale, home_dir=None):
    frames = log_frame_reader(cfg, home_dir)
    first = next(frames)

    fig, ax = plt.subplots()

    U = first['ux'].values.reshape(cfg["Ny"], cfg["Nx"])
    V = first['uy'].values.reshape(cfg["Ny"], cfg["Nx"])
    R = first['rho'].values.reshape(cfg["Ny"], cfg["Nx"])

    x, y = np.arange(cfg["Nx"]), np.arange(cfg["Ny"])
    X, Y = np.meshgrid(x, y)
    qv = ax.quiver(
            X[::qv_step, ::qv_step], Y[::qv_step, ::qv_step],
            U[::qv_step, ::qv_step], V[::qv_step, ::qv_step],
            scale       = qv_scale,
            scale_units = 'xy',
            linewidth   = 2.0,
            color       = 'black'
    )

    im = ax.imshow(R, cmap='Blues_r')
    fig.colorbar(im, ax = ax)
    ax.invert_yaxis()

    def update(frame):
        U = frame["ux"].values.reshape(cfg["Ny"], cfg["Nx"])
        V = frame["uy"].values.reshape(cfg["Ny"], cfg["Nx"])
        R = frame["rho"].values.reshape(cfg["Ny"], cfg["Nx"])

        qv.set_UVC(U[::qv_step, ::qv_step], V[::qv_step, ::qv_step])
        im.set_data(R)
        return qv, im

    a = anim.FuncAnimation(
            fig, update,
            frames   = frames,
            interval = 1000.0 / fps,
            blit     = False,
            cache_frame_data = False
    )

    writer = FFMpegWriter(fps=30, codec="h264")
    a.save(f"{cfg['name']}.mp4", writer=writer, dpi=300)
    plt.close(fig)
    display(Video(f"{cfg['name']}.mp4", embed=True))


def field_streamplot(cfg, density=0.5, home_dir=None):
    frame = next(log_frame_reader(cfg, home_dir))

    Nx = cfg["Nx"]
    Ny = cfg["Ny"]

    ux = frame["ux"].values.reshape(Ny, Nx)
    uy = frame["uy"].values.reshape(Ny, Nx)
    speed = np.sqrt(ux**2 + uy**2)

    x, y = np.meshgrid(np.arange(Nx), np.arange(Ny))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(speed, cmap="Blues_r")
    ax.streamplot(
        x, y, ux, uy,
        density=density,
        linewidth=0.5,
        color="black"
    )

    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(f"{cfg['family']}: {cfg['name']}")

    plt.tight_layout()
    plt.show()


def field_streamplot_animation(
        cfg,
        density=0.5,
        fps=30,
        home_dir=None,
        cmap="Blues_r"
):
    frames = log_frame_reader(cfg, home_dir)
    first = next(frames)

    Nx = cfg["Nx"]
    Ny = cfg["Ny"]

    x, y = np.meshgrid(np.arange(Nx), np.arange(Ny))

    ux = first["ux"].values.reshape(Ny, Nx)
    uy = first["uy"].values.reshape(Ny, Nx)
    speed = np.sqrt(ux**2 + uy**2)

    fig, ax = plt.subplots()
    im = ax.imshow(speed, cmap=cmap)

    sp = ax.streamplot(
        x, y, ux, uy,
        density=density,
        linewidth=0.5,
        color="black"
    )

    fig.colorbar(im, ax=ax)

    ax.invert_yaxis()
    ax.set_aspect("equal")


    def update(frame):
        ax.clear()
        ux = frame["ux"].values.reshape(Ny, Nx)
        uy = frame["uy"].values.reshape(Ny, Nx)
        speed = np.sqrt(ux ** 2 + uy ** 2)

        ax.imshow(speed, cmap=cmap)
        ax.streamplot(
                x, y, ux, uy,
                density=density,
                linewidth=0.5,
                color="black"
        )
        ax.invert_yaxis()
        ax.set_aspect("equal")

        return ax.collections

    ani = anim.FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=1000 / fps,
        blit=False,
        cache_frame_data=False,
    )

    writer = FFMpegWriter(
        fps=fps,
        extra_args=[
            "-crf", "18",
            "-preset", "slow",
            "-pix_fmt", "yuv420p",
        ],
    )

    filename = f"{cfg['name']}_streamplot.mp4"
    ani.save(filename, writer=writer, dpi=300)
    plt.close(fig)
    display(Video(filename, embed=True))


def profile_animation(cfg, home_dir=None):
    frames = log_frame_reader(cfg, home_dir)
    first = next(frames)

    fig, ax = plt.subplots()
    line, = ax.plot(first[cfg["disp_x"]], first[cfg["disp_y"]])
    ax.set_xlabel(cfg['disp_x'])
    ax.set_ylabel(cfg['disp_y'])

    def update(frame):
        line.set_data(frame[cfg["disp_x"]], frame[cfg["disp_y"]])
        return line,

    a = anim.FuncAnimation(
            fig, update,
            frames   = frames,
            interval = 1000.0 / fps,
            blit     = False,
            cache_frame_data = False
    )

    writer = FFMpegWriter(fps=30, codec="h264")
    a.save(f"{cfg['name']}.mp4", writer=writer, dpi=300)
    plt.close(fig)
    display(Video(f"{cfg['name']}.mp4", embed=True))


def profile_snapshots(cfg, home_dir=None):
    for frame in log_frame_reader(cfg, home_dir):
        t = frame["t"].iat[0]

        plt.plot(
                frame[cfg["disp_x"]],
                frame[cfg["disp_y"]],
                label=f"time step: {t}"
        )

    plt.xlabel(cfg["disp_x"])
    plt.ylabel(cfg["disp_y"])
    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.legend()
    plt.grid()
    plt.show()


def point_amplitude_omega_sweep(cfg, home_dir=None):
    df = log_dataframe_reader(cfg, home_dir)

    # if omega_idx is not None:
    #     values = sorted(df.omega.unique())
    #     selected = [values[i] for i in omega_idx]
    #     df = df[df.omega.isin(selected)]

    for omega, group in df.groupby(['omega']):
        plt.plot(group[cfg["disp_x"]], group[cfg["disp_y"]], label=f"omega: {omega}")
        plt.xlabel('time step')
        plt.ylabel('amplitude')

    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.legend()
    plt.grid()
    plt.show()

def point_viscosity_omega_sweep(cfg, home_dir=None):
    df = log_dataframe_reader(cfg, home_dir)

    omegas = []
    measured_nus = []
    expected_nus = []

    for omega, group in df.groupby("omega"):
        t = group.t
        A = np.abs(group.ux)
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
    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.xlabel('omega')
    plt.ylabel('nu')
    plt.legend()
    plt.grid()
    plt.show()
