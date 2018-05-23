# setup.py
from distutils.core import setup, Extension

setup(name="sample1", 
      ext_modules=[
        Extension("sample1",
                  ["vision_rlcpp02.cpp"],
                  include_dirs = ['..'],
                  language="c++",
                  )
        ]
)
