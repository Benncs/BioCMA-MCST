public_includes = include_directories(f'public/')

sources = ['src/species_descriptor.cpp']

lib_mixture = static_library('mixture', sources, include_directories: public_includes)

mixture_dependency = declare_dependency(link_with: lib_mixture, include_directories: public_includes)
if build_test
    subdir('tests')
endif