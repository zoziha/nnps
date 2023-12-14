
nnps_inc = nnps_lib.private_dir_include()
nnps_dep = declare_dependency(
    link_with: nnps_lib,
    include_directories: nnps_inc,
)


if host_machine.system() == 'windows'
    symbols_file = 'lib'+meson.project_name()+'-'+meson.project_version().split('.')[0]+'.dll.symbols'
    obj_suffix = '.obj'
else
    symbols_file = 'lib'+meson.project_name()+'.so.'+meson.project_version()+'.symbols'
    obj_suffix = '.o'
endif

nnps_obj = []
foreach source: filenames
    nnps_obj += source.replace('/', '_')+obj_suffix
endforeach

install_subdir(
    nnps_lib.full_path()+'.p',
    install_dir: 'include'/meson.project_name(),
    strip_directory: true,
    exclude_files: [
        'depscan.dd',
        meson.project_name()+'-deps.json',
        symbols_file,
        meson.project_name()+'.dat',
        nnps_obj,
    ],
)


nnps_lic = files('../LICENSE')
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

