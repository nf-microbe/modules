#!/usr/bin/env python

import argparse
import fileinput
import os
import sys


def parse_args(args=None):
    description = "Download an nf-core module into nf-microbe"
    epilog = "Example usage: python add_nfcore_module.py --module fastqc"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        help="Path to the directory containing nf-microbe/modules repository.",
    )
    parser.add_argument(
        "-m", "--module", help="Name of module to download",
    )
    return parser.parse_args(args)

# convert from modules repo to pipeline repo
def pipeline_nfcore_yml(nfcore_yml='.nf-core.yml'):
    with fileinput.FileInput(nfcore_yml, inplace=True) as file:
        for line in file:
            print(line.replace('repository_type: modules', 'repository_type: pipeline'), end='')

# convert from pipeline repo to modules repo
def modules_nfcore_yml(nfcore_yml='.nf-core.yml'):
    with fileinput.FileInput(nfcore_yml, inplace=True) as file:
        for line in file:
            print(line.replace('repository_type: pipeline', 'repository_type: modules'), end='')

# install nf-core modules
def install_nfcore_module(directory, module_name):
    os.system('nf-core modules install ' + '--dir ' + directory + " " + module_name)

# main function
def main(args=None):
    args = parse_args(args)

    pipeline_nfcore_yml(args.directory + '/.nf-core.yml')

    install_nfcore_module(
        args.directory,
        args.module
    )

    modules_nfcore_yml(args.directory + '/.nf-core.yml')

if __name__ == "__main__":
    sys.exit(main())
