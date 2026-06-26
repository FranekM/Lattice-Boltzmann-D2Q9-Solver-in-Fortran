from plots import *


swd_show = False
dwd_show = False
hbc_show = False
pf_show = False
cf_show = True


if swd_show: # shear wave decay
    configs = init("shear-wave-decay")

    config = configs['swd-field-animation']
    field_animation(config, 2, 0.005)

    # config = configs['swd-field-animation']
    # fld_streamplot(config, 1)

    # config = configs['swd-profile-animation']
    # profile_animation(config)

    config = configs['swd-profile-plots']
    profile_plot(config)

    config = configs['swd-amplitude-sweep']
    amp_omega_plot(config)

    config = configs['swd-viscosity-sweep']
    vis_omega_plot(config)


if dwd_show: # density wave decay
    configs = init("density-wave-decay")

    config = configs['dwd-field-animation']
    field_animation(config, 2, 0.005)

    config = configs['dwd-profile-plots']
    profile_plot(config)

    config = configs['dwd-amplitude-sweep']
    amp_omega_plot(config, [0, 3, 6])


if pf_show: # poiseuille flow
    configs = init('poiseuille-flow')

    config = configs['pf-field-animation']
    field_animation(config, 4, 0.005)

    config = configs['pf-profile']
    profile_plot(config)


if cf_show: # couette flow
    configs = init('couette-flow')

    config = configs['cf-field-animation']
    field_animation(config, 4, 5.0)

    config = configs['cf-profile']
    profile_plot(config)


if hbc_show: # closed box pulse
    configs = init('hard-bc')

    config = configs['box-density-pulse-animation']
    field_animation(config, 8, 0.005)

    config = configs['tube-density-pulse-animation']
    field_animation(config, 8, 0.005)



