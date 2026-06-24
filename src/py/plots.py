import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import f90nml
import math

fps = 24


def init(config_file):
    configs_raw = f90nml.read(f"../../config/{config_file}.in")['experiment']
    config = {}
    for c in configs_raw:
        c = c['cfg']
        config[c['name']] = c
    return config


def field_reader(cfg):
    reader = pd.read_csv(
            '../../' + cfg['log']["file"],
            sep       = '\t',
            names     = ['t', 'y', 'x', 'rho', 'uy', 'ux'],
            chunksize = cfg["Nx"] * cfg["Ny"]
    )

    for frame in reader:
        yield frame

def profile_reader(cfg):
    frames = pd.read_csv(
            '../../' + cfg["log"]["file"],
            sep         = '\t',
            names       = ['t', 'x', 'rho', 'uy', 'ux'],
            chunksize   = cfg["Nx"]
    )
    for frame in frames:
        yield frame


def profile_animation(cfg):
    frames = profile_reader(cfg)
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
            blit     = False
    )
    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.grid()
    plt.show()

def field_animation(cfg, qv_step, qv_scale):
    frames = field_reader(cfg)
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
            blit     = False
    )

    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.show()

def anim_fld_particles(cfg):
    # TODO correct non phisical behaviour of particles
    frames = field_reader(cfg)
    first = next(frames)

    fig, ax = plt.subplots()
    rho = first['rho'].values.reshape(cfg["Ny"], cfg["Nx"]).T

    rho_norm = rho / rho.mean()
    particle_count = (rho_norm * 5).astype(int)
    px = []
    py = []

    for y in range(cfg['Ny']):
        for x in range(cfg['Nx']):
            n = particle_count[y, x]
            px.extend(x + np.random.random(n))
            py.extend(y + np.random.random(n))

    px = np.array(px)
    py = np.array(py)
    particles = ax.scatter(
            x = px,
            y = py,
            s = 0.5,
            color = 'black'
    )

    im = ax.imshow(rho, cmap='Blues_r')
    fig.colorbar(im, ax=ax)
    ax.invert_yaxis()


    def update(frame):
        nonlocal px, py

        ux = frame["ux"].values.reshape(cfg["Ny"], cfg["Nx"]).T
        uy = frame["uy"].values.reshape(cfg["Ny"], cfg["Nx"]).T
        rho = frame["rho"].values.reshape(cfg["Ny"], cfg["Nx"]).T

        ix = np.clip(px.astype(int),0, cfg["Nx"]-1)
        iy = np.clip(py.astype(int),0, cfg["Ny"]-1)
        px += ux[iy, ix] * 50.0
        py += uy[iy, ix] * 50.0
        px %= cfg["Nx"]
        py %= cfg["Ny"]

        particles.set_offsets( np.column_stack([px, py]) )
        im.set_data(rho)
        return particles, im


    a = anim.FuncAnimation(
            fig, update,
            frames=frames,
            interval=1000.0 / fps,
            blit=False
    )
    plt.show()


def profile_plot(cfg):
    df = pd.read_csv(
            '../../' + cfg["log"]["file"],
            sep         = '\t',
            names       = ['t', 'x', 'rho', 'uy', 'ux']
    )

    for t, group in df.groupby(['t']):
        plt.plot(group[cfg["disp_x"]], group[cfg["disp_y"]], label=f"time step: {t}")
        plt.xlabel(cfg["disp_x"])
        plt.ylabel(cfg["disp_y"])

    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.legend()
    plt.grid()
    plt.show()

def amp_omega_plot(cfg, omega_idx=None):
    df = pd.read_csv(
            '../../' + cfg["log"]["file"],
            sep         = '\t',
            names       = ['omega', 't', 'rho', 'uy', 'ux']
    )

    if omega_idx is not None:
        values = sorted(df.omega.unique())
        selected = [values[i] for i in omega_idx]
        df = df[df.omega.isin(selected)]

    for omega, group in df.groupby(['omega']):
        plt.plot(group[cfg["disp_x"]], group[cfg["disp_y"]], label=f"omega: {omega}")
        plt.xlabel('time step')
        plt.ylabel('amplitude')

    plt.title(f"{cfg['family']}: {cfg['name']}")
    plt.legend()
    plt.grid()
    plt.show()

def fld_streamplot(cfg, Np=1, t_list=None):
    Nt = cfg['Nt']
    Nx = cfg['Nx']
    Ny = cfg['Ny']
    frames = field_reader(cfg)

    if t_list is not None:
        Np = len(t_list)
        selected = np.array(t_list, dtype=int)
    else:
        selected = np.linspace(0, Nt - 1, Np, dtype=int)

    cols = math.ceil(math.sqrt(Np))
    rows = math.ceil(Np / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 3*rows))
    axes = np.array(axes).flatten()

    x, y = np.meshgrid(np.arange(Nx), np.arange(Ny))

    last_frame = -1
    for i, frame_id in enumerate(selected):
        for _ in range(frame_id - last_frame - 1):
            next(frames)

        frame = next(frames)
        last_frame = frame_id
        ax = axes[i]

        ux = frame['ux'].values.reshape(Ny, Nx).T
        uy = frame['uy'].values.reshape(Ny, Nx).T
        speed = np.sqrt(ux**2 + uy**2)

        ax.imshow(speed, cmap="Blues_r")
        ax.streamplot(
                x, y, ux, uy,
                density = 0.5,
                linewidth= 0.5,
                color = 'black'
        )
        ax.set_title(f"t={frame_id}")
        ax.set_aspect("equal")

    for ax in axes[Np:]:
        ax.remove()

    plt.suptitle(f"{cfg['family']}: {cfg['name']}")
    plt.tight_layout()
    plt.show()


def vis_omega_plot(cfg):
    df = pd.read_csv(
            '../../' + cfg["log"]["file"],
            sep='\t',
            names=['omega', 't', 'rho', 'uy', 'ux']
    )

    omegas = []
    measured_nus = []
    expected_nus = []

    for omega, group in df.groupby("omega"):
        print(omega)
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
