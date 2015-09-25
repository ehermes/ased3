from numpy.distutils.core import Extension

d3fort = Extension(name = 'd3_fort',
                   sources = ['d3_fort.pyf', 'd3_fort.f90'])

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = 'ased3',
          description = 'Grimme D3 dispersion calculator for ASE',
          author = 'Eric Hermes',
          author_email = 'ehermes@wisc.edu',
          ext_modules = [d3fort])
