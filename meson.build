project('dahu',
  license: 'MIT',
  meson_version: '>= 0.60',
  version: run_command('version.py',
  			           check:true).stdout().strip(),
)


# Seek the backend
if meson.backend() != 'ninja'
  error('Ninja backend required')
endif

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

py.install_sources([
    'version.py',
],
  pure: false,    # Will be installed next to binaries
  subdir: 'dahu'  # Folder relative to site-packages to install to
)

subdir('src/dahu')
subdir('plugins')
