
def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration('ased3', parent_package, top_path)

    config.add_extension(name='_d3',
            sources=['_d3.pyf', '_d3params.f90', '_d3.f90'],
            libraries=['iomp5'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
