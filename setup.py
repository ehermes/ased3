from numpy.distutils.core import Extension

def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)

    config.add_subpackage('ased3')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(name='ased3',
          description='Grimme D3 dispersion calculator for ASE',
          author='Eric Hermes',
          author_email='ehermes@wisc.edu',
          configuration=configuration)
