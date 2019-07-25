from setuptools import setup, find_packages
import subprocess
import sys

PACKAGE_NAME = 'VarBen'
VERSION = '1.0'
DESCRIPTION = 'VarBen package'
AUTHOR = 'Fang ShuangSang & Li ZiYang'
LICENSE = 'Copyright 2018 .'


def check_java():
    p = subprocess.Popen(['java', '-version'], stderr=subprocess.PIPE)
    for line in p.stderr:
        if line.startswith('java version'):
            return True

    return False


def check_bwa():
    p = subprocess.Popen(['bwa'], stderr=subprocess.PIPE)
    for line in p.stderr:
        if line.startswith('Version:'):
            major, minor, sub = line.strip().split()[1].split('.')
            sub = sub.split('-')[0]
            sub_digit = ''.join([i for i in sub if i.isdigit()])
            if int(major) >= 0 and int(minor) >= 7 and int(sub_digit) >= 12:
                return True
    return False


def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_python():
    return sys.hexversion >= 0x20702f0


if __name__ == '__main__':
    if not check_python():
        sys.exit('Dependency problem: python >= 2.7.2 is required')
    if not check_bwa():
        sys.exit('Dependency problem: bwa >= 0.7.12 not found')
    if not check_samtools():
        sys.exit('Dependency problem: samtools >= 1.2 not found')

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    license=LICENSE,
    scripts=['bin/muteditor.py',
             'bin/sveditor.py',
             ],
    install_requires=[
        'pysam >=0.9.1.4',
        'numpy'
    ],
    packages=find_packages(),
    include_package_data=True,
)
