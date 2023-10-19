project (
    'nnps',
    'fortran',
    version: '1.3.202310',
    license: 'MIT',
    meson_version: '>=0.60.0',
    default_options : [
        'buildtype=debugoptimized',
        'fortran_std=f2008',
        'default_library=both',
    ],
)

subdir('src')
subdir('example')

nnps_lic = files('LICENSE')
install_data(
    nnps_lic,
    install_dir: join_paths(get_option('prefix'), 'share', 'licenses', meson.project_name()),
)

pkg = import('pkgconfig')
pkg.generate(
    name : meson.project_name(),
    description : 'NNPS scheme practice in CFD particle method.',
    version : meson.project_version(),
    libraries : nnps_lib,
    subdirs : meson.project_name(),
)
