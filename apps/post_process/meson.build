
src_dir = '.'
py = import('python').find_installation(pure: false,required: true)

# python_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.py ').stdout().splitlines()

python_sources = files(
    './__init__.py',
    './__main__.py',
    './io.py',
    './read_results.py',
    './properties.py',
    './initialiser.py',
)



py.install_sources(python_sources, subdir: 'biomc')