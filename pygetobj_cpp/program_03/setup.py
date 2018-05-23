# setup.py
from distutils.core import setup, Extension

setup(name="sample1", 
      ext_modules=[
        Extension("sample1",
                  ["vision_rlcpp03.cpp"],
                  include_dirs = ['..'],
                  language="c++",
                  )
        ]
)
