filenames = [
    'nnps_module.f90',
    'direct/nnps_direct1d_module.f90',
    'direct/nnps_direct2d_module.F90',
    'direct/nnps_direct3d_module.F90',
    'grid/nnps_grid1d_module.f90',
    'grid/nnps_grid2d_module.F90',
    'grid/nnps_grid3d_module.F90',
    'tree/nnps_tree1d_module.f90',
    'tree/nnps_tree1d_shape.f90',
    'tree/nnps_tree1d_binarytree.f90',
    'tree/nnps_tree2d_module.F90',
    'tree/nnps_tree2d_shape.f90',
    'tree/nnps_tree2d_quadtree.f90',
    'tree/nnps_tree3d_module.F90',
    'tree/nnps_tree3d_shape.f90',
    'tree/nnps_tree3d_octree.f90',
    'utils/nnps_spatial_hashing.F90',
    'utils/nnps_key_value.F90',
    'utils/nnps_kinds.F90',
    'utils/nnps_math.f90',
    'utils/nnps_vector.f90',
    'utils/nnps_int_vector.f90',
]

nnps_lib = library(
    meson.project_name(),
    sources: files(filenames),
    version: meson.project_version(),
    install: true,
    dependencies: dependency('openmp'),
)

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
