from distutils.core import setup, Extension

ext_modules = [Extension("pso",sources = ['src/pso.c'])]
setup(name = "pso", author = "Luckyson Khaidem",author_email = "khaidem90@gmail.com",description="Particle Swarm Optimizaiton",version = "0.1",ext_modules = ext_modules)
