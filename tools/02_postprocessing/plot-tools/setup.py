'''
Date: 2021-08-21 13:28:51
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:08:06
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(
        name='abacus_plot',
        version='1.2.0',
        packages=find_packages(),
        description='Ploting tools for ABACUS',
        author='jiyuyang',
        author_email='jiyuyang@mail.ustc.edu.cn',
        url='None',
        entry_points={'console_scripts': ['abacus-plot=abacus_plot.main:main']},
        install_requires=['matplotlib', 'numpy', 'setuptools', 'lxml'],
        python_requires='>=3.4'
    )
