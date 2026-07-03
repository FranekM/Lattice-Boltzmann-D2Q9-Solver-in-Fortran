from plots import *


swd_show = False
dwd_show = False
gwd_show = False
pf_show = False
cf_show = False
mw_show = False
sl_show = False
sl300_show = True


if swd_show: # shear wave decay
    configs = init("shear-wave-decay")

    config = configs['swd-field-animation']
    field_animation(config, 2, 0.005)

    config = configs['swd-profile-animation']
    profile_animation(config)

    config = configs['swd-profile-snapshots']
    profile_snapshots(config)

    config = configs['swd-2n-profile-snapshots']
    profile_snapshots(config)

    config = configs['swd-amplitude-sweep']
    point_amplitude_omega_sweep(config)

    config = configs['swd-viscosity-sweep']
    point_viscosity_omega_sweep(config)


if dwd_show: # density wave decay
    configs = init("density-wave-decay")

    config = configs['dwd-field-animation']
    field_animation(config, 2, 0.005)

    config = configs['dwd-profile-snapshots']
    profile_snapshots(config)

    config = configs['dwd-amplitude-sweep']
    point_amplitude_omega_sweep(config)


if gwd_show: # closed box pulse
    configs = init('gaussian-wave-decay')

    config = configs['box-density-pulse-animation']
    field_animation(config, 8, 0.005)

    config = configs['tube-density-pulse-animation']
    field_animation(config, 8, 0.005)


if pf_show: # poiseuille flow
    configs = init('poiseuille-flow')

    config = configs['pf-field-animation']
    field_animation(config, 4, 0.005)

    config = configs['pf-profile-snapshots']
    profile_snapshots(config)


if cf_show: # couette flow
    configs = init('couette-flow')

    config = configs['cf-field-animation']
    field_animation(config, 4, 0.01)

    config = configs['cf-profile-snapshots']
    profile_snapshots(config)


if mw_show: # sliding lid
    configs = init('moving-walls')

    config = configs['sliding-lid-animation']
    field_streamplot(config, t_list=[5], density=1.5)

    # probably config issues
    # config = configs['big-vortex-animation']
    # field_streamplot(config, t_list=[5], density=1.5)


if sl_show:
    configs = init('sliding-lid')

    config = configs['sliding-lid']
    field_streamplot(config, t_list=[1], density=1.5)


if sl300_show:
    configs = init('sliding-lid-300')

    config = configs['sliding-lid-300']
    field_streamplot(config, t_list=[1], density=1.5)